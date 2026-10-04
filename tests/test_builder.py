import weakref
from collections import deque
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pytest
from pysam import AlignedSegment
from pysam import AlignmentFile
from pysam import AlignmentHeader

from streampile import Pileup
from streampile import PileupRead
from streampile import PileupReadType
from streampile import StreamingPileupBuilder

from .records import HEADER
from .records import record
from .records import unmapped
from .records import write_bam

BASE = PileupReadType.base
DELETION = PileupReadType.deletion
INSERTION = PileupReadType.insertion


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


def test_builder_refuses_records_not_declared_coordinate_sorted() -> None:
    header = AlignmentHeader.from_text("@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:100\n")
    with pytest.raises(ValueError, match="Records must be coordinate sorted."):
        StreamingPileupBuilder([record("r", 10, "4M", "ACGT", header=header)])


def test_builder_refuses_records_out_of_order() -> None:
    reads = [record("b", 20, "4M", "ACGT"), record("a", 10, "4M", "ACGT")]
    with (
        StreamingPileupBuilder(reads) as builder,
        pytest.raises(ValueError, match="out of coordinate order at a"),
    ):
        builder.pileup("chr1", 30)


def test_builder_is_forward_only() -> None:
    builder = StreamingPileupBuilder([record("r", 100, "4M", "ACGT")])
    builder.pileup("chr2", 50)
    with pytest.raises(ValueError, match="Attempted to advance to chr1:100 from chr2:50."):
        builder.pileup("chr1", 100)
    with pytest.raises(ValueError, match="Attempted to advance to chr2:20 from chr2:50."):
        builder.pileup("chr2", 20)


def test_builder_refuses_unknown_contigs_and_negative_positions() -> None:
    builder = StreamingPileupBuilder([record("r", 100, "4M", "ACGT")])
    with pytest.raises(ValueError, match="Contig chr3 is not in the header."):
        builder.pileup("chr3", 1)
    with pytest.raises(ValueError, match="Position must be non-negative"):
        builder.pileup("chr1", -1)
    with pytest.raises(ValueError, match="End 1 is before start 2."):
        list(builder.columns("chr1", 2, 1))


def test_builder_refuses_pileups_once_closed() -> None:
    with StreamingPileupBuilder([record("r", 10, "4M", "ACGT")]) as builder:
        builder.pileup("chr1", 11)
    for pos in (11, 12):
        with pytest.raises(ValueError, match="The builder is closed."):
            builder.pileup("chr1", pos)
    with pytest.raises(ValueError, match="The builder is closed."):
        list(builder.columns("chr1", 12, 14))


def test_builder_returns_the_same_pileup_for_a_repeated_position() -> None:
    read = record("r", 100, "4M", "ACGT")
    builder = StreamingPileupBuilder([read])
    pileup = builder.pileup("chr1", 101)
    assert builder.pileup("chr1", 101) is pileup
    assert builder.previous_pileup is pileup
    assert pileup == Pileup("chr1", 101, [PileupRead(read, 1, 1, BASE)])


def test_builder_with_no_records_or_header_is_still_forward_only() -> None:
    with StreamingPileupBuilder([]) as builder:
        assert builder.header is None
        assert builder.pileup("chr2", 50) == Pileup("chr2", 50, [])
        assert builder.pileup("chr1", 10) == Pileup("chr1", 10, [])
        with pytest.raises(ValueError, match="Attempted to advance to chr1:5 from chr1:10."):
            builder.pileup("chr1", 5)
        with pytest.raises(ValueError, match="Attempted to advance to chr2:60 from chr1:10."):
            builder.pileup("chr2", 60)


def test_builder_checks_the_header_of_an_empty_alignment_file(tmp_path: Path) -> None:
    path = write_bam(tmp_path / "empty.bam", [])
    with AlignmentFile(str(path)) as reads, StreamingPileupBuilder(reads) as builder:
        assert builder.header is reads.header
        assert builder.pileup("chr2", 5) == Pileup("chr2", 5, [])
        with pytest.raises(ValueError, match="Contig chr3 is not in the header."):
            builder.pileup("chr3", 1)
        with pytest.raises(ValueError, match="Attempted to advance to chr1:5 from chr2:5."):
            builder.pileup("chr1", 5)
    header = AlignmentHeader.from_text("@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:100\n")
    path = write_bam(tmp_path / "queryname.bam", [], header=header)
    with (
        AlignmentFile(str(path)) as reads,
        pytest.raises(ValueError, match="Records must be coordinate sorted."),
    ):
        StreamingPileupBuilder(reads)


def test_builder_piles_up_a_deletion() -> None:
    builder = StreamingPileupBuilder([record("r", 10, "2M2D2M", "ACGT", quals=[30, 31, 32, 33])])
    assert entries(builder.pileup("chr1", 11)) == [("r", "base", 1, 1, None)]
    pileup = builder.pileup("chr1", 12)
    assert entries(pileup) == [("r", "deletion", None, 2, None)]
    assert (
        pileup.pileups[0].is_del and pileup.pileups[0].qual == 32 and pileup.pileups[0].base is None
    )
    assert entries(builder.pileup("chr1", 13)) == [("r", "deletion", None, 2, None)]
    assert entries(builder.pileup("chr1", 14)) == [("r", "base", 2, 2, None)]


def test_builder_piles_up_insertions_at_read_starts_and_ends() -> None:
    reads = [
        record("opens", 10, "2I4M", "TTACGT", quals=[10, 11, 30, 30, 30, 30]),
        record("closes", 10, "4M2I", "ACGTCC"),
        record("clipped", 10, "1S1I4M", "GTACGT"),
    ]
    builder = StreamingPileupBuilder(reads)
    at_nine = builder.pileup("chr1", 9)
    assert entries(at_nine) == [
        ("opens", "insertion", None, None, "TT"),
        ("clipped", "insertion", None, None, "T"),
    ]
    assert at_nine.pileups[0].is_ins and at_nine.pileups[0].inserted_qualities == [10, 11]
    assert at_nine.unfiltered_depth == 0
    assert entries(builder.pileup("chr1", 10)) == [
        ("opens", "base", 2, 2, None),
        ("closes", "base", 0, 0, None),
        ("clipped", "base", 2, 2, None),
    ]
    assert entries(builder.pileup("chr1", 13)) == [
        ("opens", "base", 5, 5, None),
        ("closes", "base", 3, 3, None),
        ("closes", "insertion", None, None, "CC"),
        ("clipped", "base", 5, 5, None),
    ]
    assert entries(builder.pileup("chr1", 14)) == []


def test_reads_with_no_stored_qualities_pass_every_floor() -> None:
    builder = StreamingPileupBuilder(
        [record("r", 10, "2M1D1M2I", "ACGTT", quals="*")], min_base_quality=60
    )
    at_base = builder.pileup("chr1", 11)
    assert at_base.get_query_sequences == ["C"]
    assert (at_base.get_query_qualities, at_base.filtered_depth) == ([255], 1)
    at_deletion = builder.pileup("chr1", 12)
    assert (at_deletion.pileups[0].qual, at_deletion.filtered_depth) == (255, 1)
    at_insertion = builder.pileup("chr1", 13)
    assert entries(at_insertion) == [
        ("r", "base", 2, 2, None),
        ("r", "insertion", None, None, "TT"),
    ]
    assert at_insertion.pileups[1].inserted_qualities == [255, 255]


def test_reads_with_no_stored_bases_hold_no_base_or_quality() -> None:
    pileup = StreamingPileupBuilder([record("r", 10, "2M1I1M", "*")]).pileup("chr1", 11)
    assert entries(pileup) == [("r", "base", 1, 1, None), ("r", "insertion", None, None, None)]
    assert [(entry.base, entry.qual) for entry in pileup.pileups] == [(None, None), (None, None)]
    assert pileup.pileups[1].inserted_qualities is None
    assert (pileup.unfiltered_depth, pileup.filtered_depth) == (1, 0)
    assert (pileup.get_query_sequences, pileup.get_query_qualities) == ([], [])


def test_builder_leaves_out_reads_with_no_reference_consuming_operator() -> None:
    reads = [record("clipped", 10, "2S2I", "ACGT"), record("inserted", 10, "4I", "ACGT")]
    evicted: list[AlignedSegment] = []
    with StreamingPileupBuilder(reads, tap=evicted.append) as builder:
        assert [len(pileup.pileups) for pileup in builder.columns("chr1", 8, 12)] == [0, 0, 0, 0]
    assert evicted == reads
    assert Pileup.from_alignments(reads, "chr1", 9).pileups == []


def test_builder_skips_soft_and_hard_clips() -> None:
    builder = StreamingPileupBuilder([record("r", 10, "5H2S3M1S", "TTACGA")])
    assert entries(builder.pileup("chr1", 9)) == []
    assert entries(builder.pileup("chr1", 10)) == [("r", "base", 2, 2, None)]
    assert builder.pileup("chr1", 10).get_query_sequences == ["A"]
    assert entries(builder.pileup("chr1", 12)) == [("r", "base", 4, 4, None)]
    assert entries(builder.pileup("chr1", 13)) == []


def test_builder_piles_up_reference_skips() -> None:
    builder = StreamingPileupBuilder([record("r", 10, "2M3N2M", "ACGT")], min_base_quality=0)
    columns = list(builder.columns("chr1", 10, 18))
    assert [entries(pileup) for pileup in columns] == [
        [("r", "base", 0, 0, None)],
        [("r", "base", 1, 1, None)],
        [("r", "skip", None, None, None)],
        [("r", "skip", None, None, None)],
        [("r", "skip", None, None, None)],
        [("r", "base", 2, 2, None)],
        [("r", "base", 3, 3, None)],
        [],
    ]
    assert [pileup.unfiltered_depth for pileup in columns] == [1, 1, 1, 1, 1, 1, 1, 0]
    assert [pileup.filtered_depth for pileup in columns] == [1, 1, 0, 0, 0, 1, 1, 0]
    skip = columns[2].pileups[0]
    assert skip.is_refskip and not skip.is_del and not skip.is_ins
    assert (skip.base, skip.qual) == (None, None)
    assert (columns[2].get_query_sequences, columns[2].get_query_qualities) == ([], [])


def test_builder_piles_up_an_insertion_after_a_reference_skip() -> None:
    builder = StreamingPileupBuilder([record("r", 10, "2M2N1I2M", "ACTGT")])
    assert entries(builder.pileup("chr1", 13)) == [
        ("r", "skip", None, None, None),
        ("r", "insertion", None, None, "T"),
    ]
    assert entries(builder.pileup("chr1", 14)) == [("r", "base", 3, 3, None)]


def test_builder_piles_up_both_overlapping_mates() -> None:
    first = record("pair", 10, "6M", "ACGTAC", flag=99)
    second = record("pair", 13, "6M", "TACGTA", flag=147)
    builder = StreamingPileupBuilder([first, second])
    pileup = builder.pileup("chr1", 14)
    assert [(entry.alignment.flag, entry.query_position) for entry in pileup.pileups] == [
        (99, 4),
        (147, 1),
    ]
    assert pileup.get_query_sequences == ["A", "A"]


def test_builder_moves_across_contigs() -> None:
    reads = [record("one", 10, "4M", "ACGT"), record("two", 10, "4M", "TTTT", contig="chr2")]
    evicted: list[AlignedSegment] = []
    builder = StreamingPileupBuilder(reads, tap=evicted.append)
    assert builder.pileup("chr1", 11).get_query_sequences == ["C"]
    assert builder.pileup("chr2", 11).get_query_sequences == ["T"]
    assert [read.query_name for read in evicted] == ["one"]
    builder.close()
    assert [read.query_name for read in evicted] == ["one", "two"]


@pytest.mark.parametrize(
    "flag,mapq,options,kept",
    [
        (0, 5, {"min_mapq": 20}, False),
        (0, 20, {"min_mapq": 20}, True),
        (256, 60, {}, False),
        (256, 60, {"include_secondary": True}, True),
        (2048, 60, {}, False),
        (2048, 60, {"include_supplementary": True}, True),
        (1024, 60, {}, False),
        (1024, 60, {"include_duplicate": True}, True),
        (512, 60, {}, False),
        (512, 60, {"include_qcfail": True}, True),
        (1, 60, {"proper_pairs_only": True}, False),
        (3, 60, {"proper_pairs_only": True}, True),
    ],
)
def test_builder_filters_reads(flag: int, mapq: int, options: dict[str, Any], kept: bool) -> None:
    read = record("r", 10, "4M", "ACGT", flag=flag, mapq=mapq)
    evicted: list[AlignedSegment] = []
    with StreamingPileupBuilder([read], tap=evicted.append, **options) as builder:
        assert builder.pileup("chr1", 10).unfiltered_depth == (1 if kept else 0)
    assert evicted == [read]


def test_builder_asks_a_read_filter_after_its_own_filters() -> None:
    reads = [
        record("q1", 100, "50M", "A" * 50),
        record("x2", 104, "50M", "A" * 50),
        record("q3", 108, "50M", "A" * 50, flag=1024),
        record("q4", 112, "50M", "A" * 50),
        unmapped("q5"),
    ]
    asked: list[str] = []

    def keep(read: AlignedSegment) -> bool:
        asked.append(read.query_name or "")
        return not asked[-1].startswith("x")

    evicted: list[AlignedSegment] = []
    with StreamingPileupBuilder(reads, read_filter=keep, tap=evicted.append) as builder:
        pileup = builder.pileup("chr1", 115)
        assert [entry.alignment.query_name for entry in pileup.pileups] == ["q1", "q4"]
    assert asked == ["q1", "x2", "q4"]
    assert evicted == reads


def test_builder_taps_every_read_once_in_input_order() -> None:
    reads = [
        record("long", 100, "50M", "A" * 50),
        record("short", 100, "40M", "A" * 40),
        record("filtered", 110, "10M", "A" * 10, flag=1024),
        record("later", 200, "50M", "A" * 50),
        record("other", 10, "4M", "ACGT", contig="chr2"),
        unmapped("unplaced"),
    ]
    evicted: list[AlignedSegment] = []
    with StreamingPileupBuilder(reads, tap=evicted.append) as builder:
        builder.pileup("chr1", 100)
        builder.pileup("chr1", 145)
        assert evicted == []
        builder.pileup("chr1", 200)
        assert [read.query_name for read in evicted] == ["long", "short", "filtered"]
    assert evicted == reads


def counted(reads: list[AlignedSegment], pulled: list[str]) -> Iterator[AlignedSegment]:
    """Yield reads, noting the name of each one as it is read."""
    for read in reads:
        pulled.append(read.query_name or "")
        yield read


@pytest.mark.parametrize("tapped", [False, True])
def test_closing_reads_the_rest_of_the_input_only_for_a_tap(tapped: bool) -> None:
    reads = [record(f"r{start}", start, "4M", "ACGT") for start in range(10, 60, 10)]
    pulled: list[str] = []
    evicted: list[AlignedSegment] = []
    tap = evicted.append if tapped else None
    with StreamingPileupBuilder(counted(reads, pulled), tap=tap) as builder:
        builder.pileup("chr1", 10)
        assert pulled == ["r10", "r20"]
    assert pulled == (["r10", "r20", "r30", "r40", "r50"] if tapped else ["r10", "r20"])
    assert evicted == (reads if tapped else [])


class Traceable(AlignedSegment):
    """A read that can be weakly referenced, to see when nothing holds it any more."""


@pytest.mark.parametrize("tapped", [False, True])
def test_builder_holds_passed_reads_only_for_a_tap(tapped: bool) -> None:
    reads = deque([
        record("long", 100, "100M", "A" * 100, kind=Traceable),
        record("short", 100, "10M", "A" * 10, kind=Traceable),
        record("filtered", 105, "10M", "A" * 10, flag=1024, kind=Traceable),
        record("later", 300, "4M", "ACGT", kind=Traceable),
    ])
    held = {read.query_name or "": weakref.ref(read) for read in reads}
    source = (reads.popleft() for _ in range(len(reads)))
    builder = StreamingPileupBuilder(source, tap=(lambda _: None) if tapped else None)
    builder.pileup("chr1", 120)
    alive = sorted(name for name, read in held.items() if read() is not None)
    assert alive == (["filtered", "later", "long", "short"] if tapped else ["later", "long"])


def test_builder_lets_reads_be_changed_before_they_are_written(tmp_path: Path) -> None:
    source = write_bam(
        tmp_path / "in.bam", [record("one", 100, "4M", "ACGT"), record("two", 200, "4M", "ACGT")]
    )
    with (
        AlignmentFile(str(source)) as reads,
        AlignmentFile(str(tmp_path / "out.bam"), "wb", template=reads) as sink,
        StreamingPileupBuilder(reads, tap=sink.write) as builder,
    ):
        for position, mapq in ((100, 7), (200, 9)):
            for entry in builder.pileup("chr1", position).pileups:
                entry.alignment.mapping_quality = mapq
    with AlignmentFile(str(tmp_path / "out.bam")) as written:
        assert [(read.query_name, read.mapping_quality) for read in written] == [
            ("one", 7),
            ("two", 9),
        ]


def test_columns_agree_with_pileups_at_every_position() -> None:
    reads = [
        record("a", 5, "3S10M2D5M", "GGG" + "ACGTACGTAC" + "TTTTT"),
        record("b", 7, "4M2I6M", "ACGTGGACGTAC"),
        record("c", 9, "2I8M", "TTACGTACGT"),
        record("d", 10, "3M4N3M", "ACGTAC"),
        record("e", 12, "8M", "ACGTACGT", flag=1024),
        record("f", 30, "5M", "ACGTA", contig="chr2"),
    ]
    swept = [
        entries(pileup)
        for contig, start, end in (("chr1", 0, 30), ("chr2", 25, 40))
        for pileup in StreamingPileupBuilder(reads).columns(contig, start, end)
    ]
    one_by_one = [
        entries(StreamingPileupBuilder(reads).pileup(contig, position))
        for contig, start, end in (("chr1", 0, 30), ("chr2", 25, 40))
        for position in range(start, end)
    ]
    from_alignments = [
        entries(
            Pileup.from_alignments(
                [read for read in reads if not read.is_duplicate], contig, position
            )
        )
        for contig, start, end in (("chr1", 0, 30), ("chr2", 25, 40))
        for position in range(start, end)
    ]
    assert swept == one_by_one == from_alignments
    assert sum(len(columns) for columns in swept) > 0


def test_builder_reads_an_alignment_file(tmp_path: Path) -> None:
    source = write_bam(
        tmp_path / "in.bam", [record("r", 10, "4M", "ACGT"), unmapped("u")], header=HEADER
    )
    with AlignmentFile(str(source)) as reads, StreamingPileupBuilder(reads) as builder:
        assert [pileup.get_query_sequences for pileup in builder.columns("chr1", 9, 15)] == [
            [],
            ["A"],
            ["C"],
            ["G"],
            ["T"],
            [],
        ]
