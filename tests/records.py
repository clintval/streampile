from array import array
from collections.abc import Iterable
from pathlib import Path
from typing import Literal

import pysam
from bedspec import Bed3
from bedspec import Territory
from pysam import AlignedSegment
from pysam import AlignmentFile
from pysam import AlignmentHeader

from streampile import Pileup

DATA = Path(__file__).parent / "data"

HEADER = AlignmentHeader.from_text(
    "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000\n"
)


def record(
    name: str,
    start: int,
    cigar: str,
    bases: str,
    *,
    contig: str = "chr1",
    flag: int = 0,
    mapq: int = 60,
    quals: list[int] | Literal["*"] | None = None,
    header: AlignmentHeader = HEADER,
    kind: type[AlignedSegment] = AlignedSegment,
) -> AlignedSegment:
    """A mapped read built field by field, with Q40 bases unless given, and `*` as in SAM."""
    read = kind(header)
    read.query_name = name
    read.flag = flag
    read.reference_name = contig
    read.reference_start = start
    read.mapping_quality = mapq
    read.cigarstring = cigar
    read.query_sequence = None if bases == "*" else bases
    if bases != "*" and quals != "*":
        read.query_qualities = array("B", [40] * len(bases) if quals is None else quals)
    return read


def unmapped(name: str, bases: str = "ACGT", header: AlignmentHeader = HEADER) -> AlignedSegment:
    """A read with no position at all."""
    read = AlignedSegment(header)
    read.query_name = name
    read.flag = 4
    read.reference_id = -1
    read.reference_start = -1
    read.query_sequence = bases
    read.query_qualities = array("B", [40] * len(bases))
    return read


def write_fasta(path: Path, contigs: dict[str, str]) -> Path:
    """Write and index a FASTA of the given contigs."""
    path.write_text("".join(f">{name}\n{bases}\n" for name, bases in contigs.items()))
    pysam.faidx(str(path))
    return path


def header_of(contigs: dict[str, str]) -> AlignmentHeader:
    """A coordinate-sorted header of the given contigs."""
    lines = "".join(f"@SQ\tSN:{name}\tLN:{len(bases)}\n" for name, bases in contigs.items())
    return AlignmentHeader.from_text(f"@HD\tVN:1.6\tSO:coordinate\n{lines}")


def territory(*spans: tuple[str, int, int]) -> Territory:
    """The territory of 0-based half-open `(contig, start, end)` spans."""
    return Territory(Bed3(contig, start=start, end=end) for contig, start, end in spans)


def write_bam(
    path: Path, reads: Iterable[AlignedSegment], header: AlignmentHeader = HEADER
) -> Path:
    """Write reads, already in coordinate order, to an indexed BAM."""
    with AlignmentFile(str(path), "wb", header=header) as sink:
        for read in reads:
            sink.write(read)
    pysam.index(str(path))
    return path


def entries(pileup: Pileup) -> list[tuple[str, str, int | None, int | None, str | None]]:
    """Each entry of a pileup as (read name, type, query position, next, inserted bases)."""
    return [
        (
            entry.alignment.query_name or "",
            entry.pileup_type.value,
            entry.query_position,
            entry.query_position_or_next,
            entry.inserted_bases,
        )
        for entry in pileup.pileups
    ]
