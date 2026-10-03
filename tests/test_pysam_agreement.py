from pathlib import Path

import pytest
from pysam import AlignmentFile

from streampile import StreamingPileupBuilder

from .records import DATA
from .records import record
from .records import write_bam

Entry = tuple[str, int, int | None, int | None, bool, int]


def ours(path: Path, contig: str, length: int, min_base_quality: int) -> list[list[Entry]]:
    """Every column, as (name, flag, position, position or next, deletion, inserted length)."""
    columns: list[list[Entry]] = []
    with (
        AlignmentFile(str(path)) as reads,
        StreamingPileupBuilder(
            reads,
            include_secondary=True,
            include_supplementary=True,
            include_duplicate=True,
            include_qcfail=True,
        ) as builder,
    ):
        for pileup in builder.columns(contig, 0, length):
            insertions = {
                (entry.alignment.query_name, entry.alignment.flag): entry.insertion_length
                for entry in pileup.pileups
                if entry.is_ins
            }
            columns.append(
                sorted(
                    (
                        entry.alignment.query_name or "",
                        entry.alignment.flag,
                        entry.query_position,
                        entry.query_position_or_next,
                        entry.is_del,
                        insertions.get((entry.alignment.query_name, entry.alignment.flag), 0),
                    )
                    for entry in pileup.pileups
                    if not entry.is_ins and (entry.qual or 0) >= min_base_quality
                )
            )
    return columns


def htslib(path: Path, contig: str, length: int, min_base_quality: int) -> list[list[Entry]]:
    """Every column from pysam's htslib pileup engine, refskips left out."""
    columns: list[list[Entry]] = [[] for _ in range(length)]
    with AlignmentFile(str(path)) as reads:
        for column in reads.pileup(
            contig,
            0,
            length,
            truncate=True,
            stepper="nofilter",
            min_base_quality=min_base_quality,
            max_depth=1_000_000,
            ignore_overlaps=False,
            ignore_orphans=False,
        ):
            columns[column.reference_pos] = sorted(
                (
                    entry.alignment.query_name or "",
                    entry.alignment.flag,
                    entry.query_position,
                    entry.query_position_or_next,
                    bool(entry.is_del),
                    max(entry.indel, 0),
                )
                for entry in column.pileups
                if not entry.is_refskip
            )
    return columns


@pytest.mark.parametrize("min_base_quality", [0, 30])
def test_columns_agree_with_htslib_on_the_fixture(min_base_quality: int) -> None:
    path = DATA / "reads.bam"
    for contig, length in (("chr1", 60), ("chr2", 40)):
        assert ours(path, contig, length, min_base_quality) == htslib(
            path, contig, length, min_base_quality
        )


def test_columns_agree_with_htslib_on_indels_and_skips(tmp_path: Path) -> None:
    reads = [
        record("a", 5, "3S10M2D5M", "GGG" + "ACGTACGTAC" + "TTTTT", quals=list(range(20, 38))),
        record("b", 7, "4M2I6M", "ACGTGGACGTAC"),
        record("c", 9, "2M1D1I5M", "ACGTACGT"),
        record("d", 10, "3M4N3M", "ACGTAC"),
        record("e", 12, "8M2I", "ACGTACGTTT"),
        record("f", 12, "2H8M3S", "ACGTACGTTTT", flag=1024),
    ]
    path = write_bam(tmp_path / "reads.bam", reads)
    for min_base_quality in (0, 30):
        assert ours(path, "chr1", 40, min_base_quality) == htslib(
            path, "chr1", 40, min_base_quality
        )
