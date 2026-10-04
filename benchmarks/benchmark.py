"""Time piling up every base of a territory with streampile and with pysam's htslib engine.

Each engine counts, at every base, the reads holding each base at the quality floor, the reads
with a deletion there, and the reads with an insertion after it, and the counts must agree.
Each engine runs in its own process, so its peak resident memory is its own.
"""

import argparse
import hashlib
import resource
import subprocess
import sys
import time
from collections import Counter
from collections.abc import Iterator
from pathlib import Path

from bedspec import Bed3
from bedspec import Bed3N
from bedspec import BedReader
from bedspec import Territory
from pysam import AlignmentFile
from pysam import FastaFile

from streampile import StreamingPileupBuilder
from streampile import tabulate

MIN_BASE_QUALITY = 30
MIN_MAPPING_QUALITY = 20
EXCLUDE_FLAGS = 0xF04


def read_territory(bed: Path) -> Territory:
    """The territory of a BED file."""
    with BedReader.from_path[Bed3N](bed) as features:
        return Territory(features)


def spans_in_header_order(bam: Path, territory: Territory) -> list[Bed3]:
    """The spans of a territory in the order of the BAM header's contigs, as `tabulate` walks."""
    with AlignmentFile(str(bam)) as reads:
        return sorted(territory, key=lambda span: (reads.get_tid(span.refname), span.start))


def streampile_counts(bam: Path, spans: list[Bed3]) -> Iterator[Counter[str]]:
    """Count each column of the territory from one forward sweep of the reads."""
    with (
        AlignmentFile(str(bam)) as reads,
        StreamingPileupBuilder(
            reads, min_mapq=MIN_MAPPING_QUALITY, include_qcfail=False
        ) as builder,
    ):
        for span in spans:
            for pileup in builder.columns(span.refname, span.start, span.end):
                counts: Counter[str] = Counter()
                for entry in pileup.pileups:
                    if entry.is_ins:
                        counts["+"] += 1
                        continue
                    quality = entry.qual
                    if quality is None or quality < MIN_BASE_QUALITY:
                        continue
                    counts["-" if entry.is_del else entry.base or "N"] += 1
                yield counts


def htslib_counts(bam: Path, spans: list[Bed3]) -> Iterator[Counter[str]]:
    """Count each column of the territory with pysam's pileup over each span.

    The "all" stepper drops reads by flag in htslib, but only the "samtools" stepper applies
    `min_mapping_quality`, so mapping quality is checked here.
    """
    with AlignmentFile(str(bam)) as reads:
        for span in spans:
            columns = reads.pileup(
                span.refname,
                span.start,
                span.end,
                truncate=True,
                stepper="all",
                flag_filter=EXCLUDE_FLAGS,
                min_mapping_quality=MIN_MAPPING_QUALITY,
                min_base_quality=0,
                max_depth=10_000_000,
                ignore_overlaps=False,
                ignore_orphans=False,
            )
            expected = span.start
            for column in columns:
                for _ in range(expected, column.reference_pos):
                    yield Counter()
                expected = column.reference_pos + 1
                counts: Counter[str] = Counter()
                for entry in column.pileups:
                    if entry.alignment.mapping_quality < MIN_MAPPING_QUALITY:
                        continue
                    if entry.indel > 0:
                        counts["+"] += 1
                    if entry.is_refskip:
                        continue
                    position = entry.query_position_or_next
                    qualities = entry.alignment.query_qualities_str
                    if qualities is None or position >= len(qualities):
                        continue
                    if ord(qualities[position]) - 33 < MIN_BASE_QUALITY:
                        continue
                    if entry.is_del or entry.query_position is None:
                        counts["-"] += 1
                    else:
                        sequence = entry.alignment.query_sequence or ""
                        counts[sequence[entry.query_position].upper()] += 1
                yield counts
            for _ in range(expected, span.end):
                yield Counter()


def rebuilt_counts(bam: Path, spans: list[Bed3], limit: int) -> Iterator[Counter[str]]:
    """Count evenly spaced columns, each rebuilt from every overlapping read's aligned pairs."""
    positions = [(span.refname, pos) for span in spans for pos in range(span.start, span.end)]
    step = max(len(positions) // limit, 1)
    with AlignmentFile(str(bam)) as reads:
        for contig, pos in positions[::step][:limit]:
            counts: Counter[str] = Counter()
            for read in reads.fetch(contig, pos, pos + 1):
                if read.flag & EXCLUDE_FLAGS or read.mapping_quality < MIN_MAPPING_QUALITY:
                    continue
                sequence = read.query_sequence or ""
                qualities = read.query_qualities_str or ""
                pairs = read.get_aligned_pairs(matches_only=True)  # pyright: ignore[reportUnknownMemberType]
                for query, ref in pairs:
                    if ref == pos and ord(qualities[query]) - 33 >= MIN_BASE_QUALITY:
                        counts[sequence[query]] += 1
            yield counts


def tabulate_rows(bam: Path, reference: Path, territory: Territory) -> Iterator[str]:
    """Tabulate the territory, as `streampile tabulate` does, without writing a table."""
    with AlignmentFile(str(bam)) as reads, FastaFile(str(reference)) as fasta:
        for base in tabulate(
            reads,
            fasta,
            territory,
            min_base_quality=MIN_BASE_QUALITY,
            min_mapping_quality=MIN_MAPPING_QUALITY,
        ):
            yield f"{base.depth},{base.ref_reads},{base.alts},{base.alt_reads}"


def peak_megabytes() -> float:
    """The peak resident memory of this process, in megabytes."""
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return peak / 1e6 if sys.platform == "darwin" else peak / 1e3


def run(engine: str, bam: Path, reference: Path, bed: Path, rebuild_columns: int) -> None:
    """Run one engine and print its time, peak memory, columns, and a digest of its counts."""
    territory = read_territory(bed)
    spans = spans_in_header_order(bam, territory)
    started = time.perf_counter()
    digest = hashlib.sha256()
    columns = 0
    if engine == "tabulate":
        for row in tabulate_rows(bam, reference, territory):
            digest.update(row.encode())
            columns += 1
    else:
        if engine == "rebuild":
            counts = rebuilt_counts(bam, spans, rebuild_columns)
        elif engine == "streampile":
            counts = streampile_counts(bam, spans)
        else:
            counts = htslib_counts(bam, spans)
        for column in counts:
            digest.update(repr(sorted(column.items())).encode())
            columns += 1
    seconds = time.perf_counter() - started
    print(f"{engine}\t{seconds:.1f}\t{peak_megabytes():.0f}\t{columns}\t{digest.hexdigest()[:12]}")


def main() -> None:
    """Run every engine in its own process and print a table of the results."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--ref", type=Path, required=True)
    parser.add_argument("--intervals", type=Path, required=True)
    parser.add_argument("--engine", choices=["streampile", "htslib", "tabulate", "rebuild"])
    parser.add_argument("--rebuild-columns", type=int, default=2000)
    args = parser.parse_args()
    if args.engine is not None:
        run(args.engine, args.bam, args.ref, args.intervals, args.rebuild_columns)
        return
    print("engine\tseconds\tpeak_mb\tcolumns\tdigest", flush=True)
    for engine in ("htslib", "streampile", "tabulate", "rebuild"):
        command = [sys.executable, __file__, "--bam", str(args.bam), "--ref", str(args.ref)]
        command += ["--intervals", str(args.intervals), "--engine", engine]
        command += ["--rebuild-columns", str(args.rebuild_columns)]
        _ = subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
