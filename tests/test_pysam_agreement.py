import random
from pathlib import Path

import pytest
from pysam import AlignedSegment
from pysam import AlignmentFile
from pysam import AlignmentHeader

from streampile import StreamingPileupBuilder

from .records import DATA
from .records import record
from .records import unmapped
from .records import write_bam

LENGTH = 200

RANDOM_HEADER = AlignmentHeader.from_text(
    f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:{LENGTH}\n"
)

TWO_CONTIGS = AlignmentHeader.from_text(
    f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:{LENGTH}\n@SQ\tSN:chr2\tLN:{LENGTH}\n"
)

Span = tuple[str, int, int]
"""A contig and a 0-based half-open span on it."""

Entry = tuple[str, int, str, int, int]
"""One entry: read name, flag, kind, query offset or -1 for none, and inserted length."""


def ours(path: Path, spans: list[Span], min_base_quality: int) -> list[list[Entry]]:
    """Every column of the spans from one builder, in the shape `htslib` gives.

    At a floor of 0 every entry is kept but leading insertions, which htslib never reports, and
    each of those must open an alignment that holds the next position. At a higher floor only
    bases and deletions at the floor are kept, since pysam's floor judges an insertion by its
    anchor base and a skip by its next base, where a skip here has no quality.
    """
    columns: list[list[Entry]] = []
    held: set[tuple[str, int, str, int]] = set()
    leading: list[tuple[str, int, str, int]] = []
    swept: set[tuple[str, int]] = set()
    with (
        AlignmentFile(str(path)) as reads,
        StreamingPileupBuilder(reads, exclude_flags=0) as builder,
    ):
        pileups = (
            pileup for contig, start, end in spans for pileup in builder.columns(contig, start, end)
        )
        for pileup in pileups:
            contig = pileup.reference_name
            swept.add((contig, pileup.reference_pos))
            anchored = {id(entry.alignment) for entry in pileup.pileups if not entry.is_ins}
            column: list[Entry] = []
            for entry in pileup.pileups:
                read = (entry.alignment.query_name or "", entry.alignment.flag)
                if entry.is_ins:
                    if id(entry.alignment) not in anchored:
                        leading.append((contig, pileup.reference_pos + 1, *read))
                    elif min_base_quality == 0:
                        offset = entry.insertion_offset
                        assert offset is not None
                        column.append((*read, "insertion", offset, entry.insertion_length))
                    continue
                held.add((contig, pileup.reference_pos, *read))
                if entry.is_refskip:
                    if min_base_quality == 0:
                        column.append((*read, "skip", -1, 0))
                elif min_base_quality == 0 or (entry.qual or 0) >= min_base_quality:
                    offset = (
                        -1 if entry.query_position_or_next is None else entry.query_position_or_next
                    )
                    column.append((*read, entry.pileup_type.value, offset, 0))
            columns.append(sorted(column))
    assert all(read in held for read in leading if read[:2] in swept)
    return columns


def htslib(path: Path, spans: list[Span], min_base_quality: int) -> list[list[Entry]]:
    """Every column of the spans from pysam's htslib pileup engine, in our conventions.

    A skip has no offset, and a deletion that no base follows has none, where htslib gives the
    query length. An insertion is its own entry, starting after the base or at the next base.
    """
    return [
        column
        for contig, start, end in spans
        for column in htslib_span(path, contig, start, end, min_base_quality)
    ]


def htslib_span(
    path: Path, contig: str, first: int, end: int, min_base_quality: int
) -> list[list[Entry]]:
    """Every column of one span from pysam's htslib pileup engine, in our conventions."""
    columns: list[list[Entry]] = [[] for _ in range(first, end)]
    with AlignmentFile(str(path)) as reads:
        for column in reads.pileup(
            contig,
            first,
            end,
            truncate=True,
            stepper="nofilter",
            min_base_quality=min_base_quality,
            max_depth=1_000_000,
            ignore_overlaps=False,
            ignore_orphans=False,
        ):
            entries: list[Entry] = []
            for entry in column.pileups:
                read = (entry.alignment.query_name or "", entry.alignment.flag)
                offset = entry.query_position_or_next
                if entry.is_refskip:
                    if min_base_quality == 0:
                        entries.append((*read, "skip", -1, 0))
                elif entry.is_del:
                    ends = offset >= (entry.alignment.infer_query_length() or 0)
                    entries.append((*read, "deletion", -1 if ends else offset, 0))
                else:
                    entries.append((*read, "base", offset, 0))
                if min_base_quality == 0 and entry.indel > 0:
                    start = offset if entry.is_del else offset + 1
                    entries.append((*read, "insertion", start, entry.indel))
            columns[column.reference_pos - first] = sorted(entries)
    return columns


@pytest.mark.parametrize("min_base_quality", [0, 30])
def test_columns_agree_with_htslib_on_the_fixture(min_base_quality: int) -> None:
    path = DATA / "reads.bam"
    spans = [("chr1", 0, 60), ("chr2", 0, 40)]
    assert ours(path, spans, min_base_quality) == htslib(path, spans, min_base_quality)


def test_columns_agree_with_htslib_on_indels_and_skips(tmp_path: Path) -> None:
    reads = [
        record("a", 5, "3S10M2D5M", "GGG" + "ACGTACGTAC" + "TTTTT", quals=list(range(20, 38))),
        record("b", 7, "4M2I6M", "ACGTGGACGTAC"),
        record("c", 9, "2M1D1I5M", "ACGTACGT"),
        record("d", 10, "3M4N3M", "ACGTAC"),
        record("e", 12, "8M2I", "ACGTACGTTT"),
        record("f", 12, "2H8M3S", "ACGTACGTTTT", flag=1024),
        record("g", 14, "2M2N1I2M", "ACTGT", quals=[40, 40, 20, 10, 40]),
        record("h", 15, "3M2D", "ACG"),
        record("i", 15, "1D1I3M", "TACG", quals=[10, 40, 40, 40]),
        record("j", 16, "1D3M", "ACG"),
        record("k", 16, "4H1I3M", "TACG"),
        record("l", 16, "2S2I", "ACGT"),
        record("m", 17, "2M1D2M1I", "ACGTA", quals="*"),
        record("n", 17, "2M1D2M1I", "*"),
    ]
    path = write_bam(tmp_path / "reads.bam", reads)
    for min_base_quality in (0, 30):
        spans = [("chr1", 0, 40)]
        assert ours(path, spans, min_base_quality) == htslib(path, spans, min_base_quality)


def random_cigar(rng: random.Random) -> list[tuple[int, str]]:
    """CIGAR operators of every kind, which may open or close with an I, D, or N, or place none."""
    if rng.random() < 0.03:
        return rng.choice([[(4, "S")], [(3, "I")], [(2, "S"), (2, "I")], [(1, "S"), (2, "I")]])
    operators: list[tuple[int, str]] = []
    if rng.random() < 0.2:
        operators.append((rng.randint(1, 3), "H"))
    if rng.random() < 0.3:
        operators.append((rng.randint(1, 3), "S"))
    if rng.random() < 0.15:
        operators.append((rng.randint(1, 3), "I"))
    if rng.random() < 0.1:
        operators.append((rng.randint(1, 2), rng.choice("DN")))
    operators.append((rng.randint(1, 6), rng.choice("M=X")))
    for _ in range(rng.randint(0, 5)):
        operator = rng.choices("MIDNP=X", weights=[5, 3, 3, 1, 1, 1, 1])[0]
        operators.append((rng.randint(1, 4), operator))
    if rng.random() < 0.8:
        operators.append((rng.randint(1, 6), "M"))
    if rng.random() < 0.1:
        operators.append((rng.randint(1, 3), rng.choice("IDN")))
    if rng.random() < 0.3:
        operators.append((rng.randint(1, 3), "S"))
    if rng.random() < 0.2:
        operators.append((rng.randint(1, 3), "H"))
    return operators


def random_reads(
    seed: int,
    count: int = 25,
    header: AlignmentHeader = RANDOM_HEADER,
    contigs: tuple[str, ...] = ("chr1",),
) -> list[AlignedSegment]:
    """Reads with random CIGARs, flags, bases, and qualities, some with no QUAL or no SEQ."""
    rng = random.Random(seed)
    reads: list[AlignedSegment] = []
    for index in range(count):
        operators = random_cigar(rng)
        length = sum(length for length, operator in operators if operator in "MIS=X")
        bases = "".join(rng.choice("ACGTN") for _ in range(length))
        missing = rng.random()
        reads.append(
            record(
                f"r{index}",
                rng.randint(0, LENGTH - 60),
                "".join(f"{length}{operator}" for length, operator in operators),
                "*" if missing < 0.03 else bases,
                flag=rng.choice([0, 4, 16, 99, 147, 256, 512, 1024, 2048]),
                mapq=rng.choice([0, 10, 60]),
                quals="*" if missing < 0.08 else [rng.choice([2, 10, 20, 30, 40]) for _ in bases],
                contig=rng.choice(contigs),
                header=header,
            )
        )
    return sorted(reads, key=lambda read: (read.reference_id, read.reference_start))


@pytest.mark.parametrize("seed", range(50))
def test_columns_agree_with_htslib_on_random_reads(tmp_path: Path, seed: int) -> None:
    path = write_bam(tmp_path / "reads.bam", random_reads(seed), header=RANDOM_HEADER)
    spans = [("chr1", 0, LENGTH)]
    for min_base_quality in (0, 13, 30):
        assert ours(path, spans, min_base_quality) == htslib(path, spans, min_base_quality)


def gapped_spans(rng: random.Random) -> list[Span]:
    """Spans of both contigs, in order, some one base long, with gaps between some of them."""
    spans: list[Span] = []
    for contig in ("chr1", "chr2"):
        position = rng.randint(0, 20)
        while position < LENGTH:
            end = min(position + rng.choice([1, 1, 3, 10, 40]), LENGTH)
            spans.append((contig, position, end))
            position = end + rng.choice([0, 1, 5, 30, 80])
    return spans


@pytest.mark.parametrize("seed", range(50))
def test_gapped_spans_of_two_contigs_agree_with_htslib(tmp_path: Path, seed: int) -> None:
    reads = random_reads(seed, count=60, header=TWO_CONTIGS, contigs=("chr1", "chr2"))
    unplaced = [unmapped(f"u{index}", header=TWO_CONTIGS) for index in range(3)]
    path = write_bam(tmp_path / "reads.bam", reads + unplaced, header=TWO_CONTIGS)
    spans = gapped_spans(random.Random(seed))
    for min_base_quality in (0, 13, 30):
        assert ours(path, spans, min_base_quality) == htslib(path, spans, min_base_quality)


def test_the_fixture_bam_holds_the_fixture_sam() -> None:
    with (
        AlignmentFile(str(DATA / "reads.sam")) as sam,
        AlignmentFile(str(DATA / "reads.bam")) as bam,
    ):
        assert [read.to_string() for read in sam] == [read.to_string() for read in bam]
