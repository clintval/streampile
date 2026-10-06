import weakref
from array import array
from collections import Counter
from collections import deque
from collections.abc import Iterator
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

import pysam
import pytest
from pysam import AlignedSegment
from pysam import AlignmentFile
from pysam import AlignmentHeader

from streampile import Pileup
from streampile import PileupRead
from streampile import PileupReadType
from streampile import StreamingPileupBuilder
from streampile._pileup import BASE

from .records import HEADER
from .records import entries
from .records import record
from .records import unmapped
from .records import write_bam


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
    assert pileup == Pileup("chr1", 101, (PileupRead(read, 1, 1, BASE),))


def test_builder_with_no_records_or_header_is_still_forward_only() -> None:
    with StreamingPileupBuilder([]) as builder:
        assert builder.header is None
        assert builder.pileup("chr2", 50) == Pileup("chr2", 50, ())
        assert builder.pileup("chr1", 10) == Pileup("chr1", 10, ())
        with pytest.raises(ValueError, match="Attempted to advance to chr1:5 from chr1:10."):
            builder.pileup("chr1", 5)
        with pytest.raises(ValueError, match="Attempted to advance to chr2:60 from chr1:10."):
            builder.pileup("chr2", 60)


def test_builder_checks_the_header_of_an_empty_alignment_file(tmp_path: Path) -> None:
    path = write_bam(tmp_path / "empty.bam", [])
    with AlignmentFile(str(path)) as reads, StreamingPileupBuilder(reads) as builder:
        assert builder.header is reads.header
        assert builder.pileup("chr2", 5) == Pileup("chr2", 5, ())
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


def test_builder_piles_up_a_deletion_no_base_follows() -> None:
    reads = [record("r", 10, "3M2D", "ACG"), record("s", 10, "4M", "ACGT")]
    builder = StreamingPileupBuilder(reads, min_base_quality=0)
    pileup = builder.pileup("chr1", 13)
    assert entries(pileup) == [("r", "deletion", None, None, None), ("s", "base", 3, 3, None)]
    assert pileup.pileups[0].qual is None
    assert (pileup.unfiltered_depth, pileup.filtered_depth) == (2, 1)
    assert entries(builder.pileup("chr1", 14)) == [("r", "deletion", None, None, None)]
    assert entries(builder.pileup("chr1", 15)) == []


@pytest.mark.parametrize(
    "cigar,bases,expected",
    [
        ("1D3M", "ACG", [("r", "deletion", None, 0, None)]),
        ("3D4M", "ACGT", [("r", "deletion", None, 0, None)]),
        (
            "1D1I3M",
            "TACG",
            [("r", "deletion", None, 0, None), ("r", "insertion", None, None, "T")],
        ),
    ],
)
def test_builder_piles_up_a_read_that_opens_with_a_deletion(
    cigar: str, bases: str, expected: list[tuple[str, str, int | None, int | None, str | None]]
) -> None:
    read = record("r", 10, cigar, bases, quals=[25] + [40] * (len(bases) - 1))
    builder = StreamingPileupBuilder([read])
    assert entries(builder.pileup("chr1", 9)) == []
    pileup = builder.pileup("chr1", 10)
    assert entries(pileup) == expected
    assert pileup.pileups[0].qual == 25


def test_builder_piles_up_an_opening_insertion_after_a_hard_clip() -> None:
    builder = StreamingPileupBuilder([record("r", 10, "4H1I3M", "TACG")])
    assert entries(builder.pileup("chr1", 9)) == [("r", "insertion", None, None, "T")]
    assert entries(builder.pileup("chr1", 10)) == [("r", "base", 1, 1, None)]


def test_builder_counts_a_column_of_every_kind_of_entry() -> None:
    reads = [
        record("base", 10, "4M", "ACGT"),
        record("lowdel", 10, "2M1D2M", "ACTA", quals=[40, 40, 5, 40]),
        record("del", 10, "2M1D2M", "ACTA"),
        record("closes", 10, "3M1I1M", "ACGTA"),
        record("opens", 13, "1I3M", "TACG"),
    ]
    pileup = StreamingPileupBuilder(reads).pileup("chr1", 12)
    assert entries(pileup) == [
        ("base", "base", 2, 2, None),
        ("lowdel", "deletion", None, 2, None),
        ("del", "deletion", None, 2, None),
        ("closes", "base", 2, 2, None),
        ("closes", "insertion", None, None, "T"),
        ("opens", "insertion", None, None, "T"),
    ]
    assert (len(pileup.pileups), pileup.unfiltered_depth, pileup.filtered_depth) == (6, 4, 3)
    assert pileup.bases == ["G", "G"]


def test_builder_floor_leaves_bases_under_it_out_of_the_views_only() -> None:
    reads = [
        record(f"q{quality}", 100, "50M", "A" * 50, quals=[quality] * 50)
        for quality in (19, 20, 21)
    ]
    pileup = StreamingPileupBuilder(reads, min_base_quality=20).pileup("chr1", 104)
    assert (pileup.filtered_depth, pileup.qualities) == (2, [20, 21])
    assert pileup.bases == ["A", "A"]
    assert [entry.alignment.query_name for entry in pileup.pileups] == ["q19", "q20", "q21"]
    assert pileup.unfiltered_depth == 3


def test_builder_piles_up_a_crowd_of_reads_at_one_start() -> None:
    reads = [
        record(f"{base}{index}", 5, "10M", base * 10)
        for base, count in (("A", 5), ("C", 4), ("G", 3), ("T", 2), ("N", 1))
        for index in range(count)
    ]
    pileup = StreamingPileupBuilder(reads).pileup("chr1", 5)
    assert Counter(pileup.bases) == {"A": 5, "C": 4, "G": 3, "T": 2, "N": 1}
    assert pileup.unfiltered_depth == 15
    assert [entry.alignment.query_name for entry in pileup.pileups] == [
        read.query_name for read in reads
    ]


def test_builder_piles_up_read_through_pairs_and_reverse_reads_as_aligned() -> None:
    reads = [
        record("pair", 99, "10M", "ACGTACGTAC", flag=83, quals=[35] * 10),
        record("pair", 100, "10M", "CGTACGTACG", flag=163, quals=[35] * 10),
    ]
    columns = StreamingPileupBuilder(reads).columns("chr1", 98, 111)
    assert [pileup.unfiltered_depth for pileup in columns] == [0, 1] + [2] * 9 + [1, 0]
    reverse = StreamingPileupBuilder(reads).pileup("chr1", 103).pileups[0]
    assert reverse.alignment.is_reverse
    assert (reverse.query_position, reverse.base, reverse.qual) == (4, "A", 35)


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
    assert [entry.insertion_offset for entry in at_nine.pileups] == [0, 1]
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
    assert at_base.bases == ["C"]
    assert (at_base.qualities, at_base.filtered_depth) == ([255], 1)
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
    assert (pileup.bases, pileup.qualities) == ([], [])


def test_builder_leaves_out_reads_with_no_reference_consuming_operator() -> None:
    reads = [record("clipped", 10, "2S2I", "ACGT"), record("inserted", 10, "4I", "ACGT")]
    evicted: list[AlignedSegment] = []
    with StreamingPileupBuilder(reads, tap=evicted.append) as builder:
        assert [len(pileup.pileups) for pileup in builder.columns("chr1", 8, 12)] == [0, 0, 0, 0]
    assert evicted == reads
    assert Pileup.from_alignments(reads, "chr1", 9).pileups == ()


def test_builder_skips_soft_and_hard_clips() -> None:
    builder = StreamingPileupBuilder([record("r", 10, "5H2S3M1S", "TTACGA")])
    assert entries(builder.pileup("chr1", 9)) == []
    assert entries(builder.pileup("chr1", 10)) == [("r", "base", 2, 2, None)]
    assert builder.pileup("chr1", 10).bases == ["A"]
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
    assert (columns[2].bases, columns[2].qualities) == ([], [])


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
    assert pileup.bases == ["A", "A"]


def test_builder_moves_across_contigs() -> None:
    reads = [record("one", 10, "4M", "ACGT"), record("two", 10, "4M", "TTTT", contig="chr2")]
    evicted: list[AlignedSegment] = []
    builder = StreamingPileupBuilder(reads, tap=evicted.append)
    assert builder.pileup("chr1", 11).bases == ["C"]
    assert builder.pileup("chr2", 11).bases == ["T"]
    assert [read.query_name for read in evicted] == ["one"]
    builder.close()
    assert [read.query_name for read in evicted] == ["one", "two"]


def test_builder_floors_bases_at_13_and_leaves_out_qc_fail_reads_by_default() -> None:
    reads = [
        record("q12", 10, "4M", "ACGT", quals=[12] * 4),
        record("q13", 10, "4M", "ACGT", quals=[13] * 4),
        record("qcfail", 10, "4M", "ACGT", flag=512),
    ]
    pileup = StreamingPileupBuilder(reads).pileup("chr1", 10)
    assert [entry.alignment.query_name for entry in pileup.pileups] == ["q12", "q13"]
    assert (pileup.filtered_depth, pileup.qualities) == (1, [13])


@pytest.mark.parametrize(
    "flag,mapq,options,kept",
    [
        (0, 5, {"min_mapping_quality": 20}, False),
        (0, 20, {"min_mapping_quality": 20}, True),
        (256, 60, {}, False),
        (256, 60, {"exclude_flags": 0xE00}, True),
        (2048, 60, {}, False),
        (2048, 60, {"exclude_flags": 0x700}, True),
        (1024, 60, {}, False),
        (1024, 60, {"exclude_flags": 0xB00}, True),
        (512, 60, {}, False),
        (512, 60, {"exclude_flags": 0xD00}, True),
        (512, 60, {"exclude_flags": 0xE00}, False),
        (16, 60, {"exclude_flags": 0x10}, False),
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


def test_builder_taps_a_read_once_it_has_passed_the_read_and_every_one_before() -> None:
    reads = [
        record("first", 10, "4M", "ACGT"),
        record("unmapped", 10, "*", "ACGT", flag=4),
        record("second", 12, "4M", "ACGT"),
    ]
    evicted: list[AlignedSegment] = []
    builder = StreamingPileupBuilder(reads, tap=evicted.append)
    pileup = builder.pileup("chr1", 13)
    assert [entry.alignment.query_name for entry in pileup.pileups] == ["first", "second"]
    assert evicted == []
    builder.pileup("chr1", 14)
    assert [read.query_name for read in evicted] == ["first", "unmapped"]


def test_builder_taps_every_read_when_its_block_raises() -> None:
    reads = [record(f"r{start}", start, "4M", "ACGT") for start in (10, 20, 30)]
    evicted: list[AlignedSegment] = []
    with pytest.raises(RuntimeError), StreamingPileupBuilder(reads, tap=evicted.append) as builder:
        builder.pileup("chr1", 10)
        raise RuntimeError
    assert evicted == reads


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
        record("f", 14, "3M2D", "ACG"),
        record("g", 30, "5M", "ACGTA", contig="chr2"),
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
        assert [pileup.bases for pileup in builder.columns("chr1", 9, 15)] == [
            [],
            ["A"],
            ["C"],
            ["G"],
            ["T"],
            [],
        ]


READ_LENGTH = 50


def pair(
    name: str, start1: int, start2: int, *, reverse1: bool = False, reverse2: bool = True
) -> list[AlignedSegment]:
    """Two 50-base mates at 0-based starts with mate fields as htsjdk sets them, sorted."""

    def five_prime(start: int, reverse: bool) -> int:
        return start + READ_LENGTH - 1 if reverse else start

    first, second = five_prime(start1, reverse1), five_prime(start2, reverse2)
    insert = second - first + (1 if second >= first else -1)
    reads: list[AlignedSegment] = []
    for start, mate_start, reverse, mate_reverse, flag, tlen, base in (
        (start1, start2, reverse1, reverse2, 65, insert, "A"),
        (start2, start1, reverse2, reverse1, 129, -insert, "C"),
    ):
        flag |= (16 if reverse else 0) | (32 if mate_reverse else 0)
        read = record(name, start, f"{READ_LENGTH}M", base * READ_LENGTH, flag=flag)
        read.next_reference_id = read.reference_id
        read.next_reference_start = mate_start
        read.template_length = tlen
        reads.append(read)
    return sorted(reads, key=lambda read: read.reference_start)


def test_builder_leaves_out_reads_on_a_previous_contig_or_ending_just_before() -> None:
    reads = [record("prev-contig", 54, "50M", "A" * 50) for _ in range(5)]
    reads += [record("abutting", 4, "50M", "A" * 50, contig="chr2") for _ in range(5)]
    reads += [
        record("here", 54, "50M", base * 50, contig="chr2")
        for base, count in (("A", 5), ("C", 4), ("G", 3), ("T", 2), ("N", 1))
        for _ in range(count)
    ]
    pileup = StreamingPileupBuilder(reads).pileup("chr2", 54)
    assert pileup.unfiltered_depth == 15
    assert {entry.alignment.query_name for entry in pileup.pileups} == {"here"}
    assert Counter(pileup.bases) == {"A": 5, "C": 4, "G": 3, "T": 2, "N": 1}


def test_builder_piles_up_every_edge_case_of_indels() -> None:
    reads = [
        record("q1", 100, "10M2D40M", "A" * 50),
        record("q2", 100, "10M2I38M", "C" * 50),
        record("q3", 100, "31M9I10M", "G" * 50),
        record("q4", 100, "30M9D20M", "T" * 50),
        record("q5", 140, "10I40M", "N" * 50),
        record("q6", 200, "47M3S", "N" * 50),
    ]
    builder = StreamingPileupBuilder(reads)

    def named(pileup: Pileup, kind: PileupReadType) -> list[str | None]:
        return [entry.alignment.query_name for entry in pileup.pileups if entry.pileup_type is kind]

    assert len(builder.pileup("chr1", 104).bases) == 4
    before = builder.pileup("chr1", 109)
    assert (before.unfiltered_depth, len(before.pileups), sorted(before.bases)) == (
        4,
        5,
        list("ACGT"),
    )
    assert named(before, PileupReadType.insertion) == ["q2"]
    deleted = builder.pileup("chr1", 110)
    assert (deleted.unfiltered_depth, len(deleted.pileups), sorted(deleted.bases)) == (
        4,
        4,
        list("CGT"),
    )
    assert named(deleted, PileupReadType.deletion) == ["q1"]
    bigger = builder.pileup("chr1", 130)
    assert (bigger.unfiltered_depth, len(bigger.pileups), sorted(bigger.bases)) == (
        4,
        5,
        list("ACG"),
    )
    assert named(bigger, PileupReadType.insertion) == ["q3"]
    assert named(bigger, PileupReadType.deletion) == ["q4"]
    leading = builder.pileup("chr1", 139)
    assert (leading.unfiltered_depth, len(leading.pileups), sorted(leading.bases)) == (
        4,
        5,
        list("ACGT"),
    )
    assert named(leading, PileupReadType.insertion) == ["q5"]
    assert leading.pileups[4].insertion_offset == 0
    clipped = builder.pileup("chr1", 246)
    assert (clipped.unfiltered_depth, len(clipped.pileups), clipped.bases) == (1, 1, ["N"])
    past = builder.pileup("chr1", 247)
    assert (past.unfiltered_depth, len(past.pileups)) == (0, 0)


def test_builder_piles_up_only_reads_of_mapped_pairs_with_a_read_filter() -> None:
    half_mapped = pair("q2", 100, 100, reverse2=False)
    half_mapped[0].flag = 1 | 8 | 64
    half_mapped[1].flag = 1 | 4 | 128
    reads = [record("q1", 100, "50M", "A" * 50), *half_mapped, *pair("q3", 100, 299)]

    def mapped_pair(read: AlignedSegment) -> bool:
        return read.is_paired and not read.is_unmapped and not read.mate_is_unmapped

    pileup = StreamingPileupBuilder(reads, read_filter=mapped_pair).pileup("chr1", 104)
    assert pileup.unfiltered_depth == 1
    kept = pileup.pileups[0].alignment
    assert (kept.query_name, kept.is_read1) == ("q3", True)


def test_builder_keeps_a_fragment_and_positions_outside_the_insert_of_an_fr_pair() -> None:
    fragment = StreamingPileupBuilder([record("q1", 99, "50M", "A" * 50)])
    assert [fragment.pileup("chr1", pos).unfiltered_depth for pos in (99, 148)] == [1, 1]
    builder = StreamingPileupBuilder(pair("q2", 100, 99))
    assert [builder.pileup("chr1", pos).unfiltered_depth for pos in (99, 100, 148, 149)] == [
        1,
        2,
        2,
        1,
    ]


def test_builder_keeps_every_position_of_a_pair_with_the_reverse_read_starting_later() -> None:
    builder = StreamingPileupBuilder(pair("q2", 100, 99, reverse1=True, reverse2=False))
    assert [builder.pileup("chr1", pos).unfiltered_depth for pos in (99, 100, 148, 149)] == [
        1,
        2,
        2,
        1,
    ]


def test_builder_composes_a_read_filter_with_an_entry_filter() -> None:
    reads = [
        record(name, start, "50M", "A" * 50)
        for name, start in (("q1", 100), ("x2", 104), ("q3", 108), ("q4", 112))
    ]

    def keep(read: AlignedSegment) -> bool:
        return not (read.query_name or "").startswith("x")

    pileup = StreamingPileupBuilder(reads, read_filter=keep).pileup("chr1", 114)
    kept = [
        entry.alignment.query_name
        for entry in pileup.pileups
        if entry.query_position is not None and entry.query_position > 5
    ]
    assert kept == ["q1", "q3"]


def test_entries_report_offsets_in_alignment_order_on_both_strands() -> None:
    reads = pair("q1", 100, 200)
    for read in reads:
        read.query_qualities = array("B", [35] * READ_LENGTH)
    builder = StreamingPileupBuilder(reads)
    seen: list[tuple[str | None, int | None, int | None, bool]] = []
    for pos in (104, 204):
        pileup = builder.pileup("chr1", pos)
        assert pileup.unfiltered_depth == 1
        entry = pileup.pileups[0]
        seen.append((entry.base, entry.qual, entry.query_position, entry.alignment.is_reverse))
    assert seen == [("A", 35, 4, False), ("C", 35, 4, True)]


def test_entries_measure_their_distances_to_both_fragment_ends() -> None:
    reads = pair("q1", 100, 150)
    for read in reads:
        read.set_tag("MC", f"{READ_LENGTH}M")  # pyright: ignore[reportUnknownMemberType]
    distances: list[tuple[bool, int | None, int | None]] = []
    with StreamingPileupBuilder(reads) as builder:
        for pos in (100, 149, 150, 199):
            for entry in builder.pileup("chr1", pos).pileups:
                distances.append((
                    entry.alignment.is_reverse,
                    entry.five_prime_distance,
                    entry.template_end_distance,
                ))
    assert distances == [
        (False, 0, 99),
        (False, 49, 50),
        (True, 49, 50),
        (True, 0, 99),
    ]


def test_a_read_needs_its_mate_cigar_for_its_template_end() -> None:
    reads = pair("q1", 100, 120)
    reads[0].template_length = 7
    pileup = StreamingPileupBuilder(reads).pileup("chr1", 130)
    missing = "Read q1 has no MC tag to find its mate's 5' end with."
    for entry in pileup.pileups:
        with pytest.raises(ValueError, match=missing):
            _ = entry.template_end_distance
    reads[0].set_tag("MC", "4Q")  # pyright: ignore[reportUnknownMemberType]
    forward = StreamingPileupBuilder(reads).pileup("chr1", 130).pileups[0]
    with pytest.raises(ValueError, match="Read q1 has an invalid MC tag: 4Q."):
        _ = forward.template_end_distance


def test_an_entry_made_by_hand_measures_its_distances_from_its_read() -> None:
    reads = pair("q1", 100, 150)
    for read in reads:
        read.set_tag("MC", f"{READ_LENGTH}M")  # pyright: ignore[reportUnknownMemberType]
    entry = PileupRead(reads[0], 4, 4, BASE)
    assert (entry.five_prime_distance, entry.template_end_distance) == (4, 95)
    reverse = PileupRead(reads[1], 4, 4, BASE)
    assert (reverse.five_prime_distance, reverse.template_end_distance) == (45, 54)
    deletion = PileupRead(reads[0], None, 4, PileupReadType.deletion)
    assert (deletion.five_prime_distance, deletion.template_end_distance) == (None, None)


def test_a_pileup_is_a_snapshot_that_outlives_the_builder_moving_on() -> None:
    reads = [record("one", 10, "4M", "ACGT"), record("two", 20, "4M", "TTTT")]
    with StreamingPileupBuilder(reads) as builder:
        first = builder.pileup("chr1", 11)
        later = builder.pileup("chr1", 21)
    assert (first.bases, [entry.alignment for entry in first.pileups]) == (["C"], [reads[0]])
    assert (later.bases, later.unfiltered_depth) == (["T"], 1)


def test_the_options_of_a_builder_are_read_only() -> None:
    def keep(_read: AlignedSegment) -> bool:
        return True

    tapped: list[AlignedSegment] = []
    builder = StreamingPileupBuilder(
        [record("r", 10, "4M", "ACGT")],
        min_mapping_quality=5,
        exclude_flags=0x400,
        min_base_quality=20,
        proper_pairs_only=True,
        read_filter=keep,
        tap=tapped.append,
    )
    assert (
        builder.min_mapping_quality,
        builder.exclude_flags,
        builder.min_base_quality,
        builder.proper_pairs_only,
        builder.read_filter,
        builder.tap,
        builder.previous_pileup,
    ) == (5, 0x400, 20, True, keep, tapped.append, None)
    with pytest.raises(AttributeError):
        builder.min_base_quality = 30  # type: ignore[misc]  # pyright: ignore[reportAttributeAccessIssue]  # ty: ignore[invalid-assignment]


def test_closing_again_after_a_tap_raises_hands_over_the_rest() -> None:
    reads = [record(f"r{start}", start, "4M", "ACGT") for start in (10, 20, 30)]
    tapped: list[str] = []

    def tap(read: AlignedSegment) -> None:
        if read.query_name == "r10" and "failed" not in tapped:
            tapped.append("failed")
            raise OSError("disk full")
        tapped.append(read.query_name or "")

    builder = StreamingPileupBuilder(reads, tap=tap)
    builder.pileup("chr1", 10)
    with pytest.raises(OSError, match="disk full"):
        builder.close()
    with pytest.raises(ValueError, match="The builder is closed."):
        builder.pileup("chr1", 40)
    builder.close()
    builder.close()
    assert tapped == ["failed", "r20", "r30"]


def test_a_failed_advance_is_retried_and_never_returns_a_stale_pileup() -> None:
    reads = [record("one", 10, "4M", "ACGT"), record("two", 20, "4M", "GGGG")]
    failures = ["disk full"]

    def tap(_read: AlignedSegment) -> None:
        if failures:
            raise OSError(failures.pop())

    builder = StreamingPileupBuilder(reads, tap=tap)
    assert builder.pileup("chr1", 10).bases == ["A"]
    with pytest.raises(OSError, match="disk full"):
        builder.pileup("chr1", 20)
    assert builder.pileup("chr1", 20).bases == ["G"]


def test_exceptions_of_the_records_read_filter_and_tap_reach_the_caller_unchanged() -> None:
    class CallbackError(Exception):
        pass

    def broken() -> Iterator[AlignedSegment]:
        yield record("one", 10, "4M", "ACGT")
        raise CallbackError("source")

    def refuse(_read: AlignedSegment) -> bool:
        raise CallbackError("filter")

    def tap(_read: AlignedSegment) -> None:
        raise CallbackError("tap")

    reads = [record("one", 10, "4M", "ACGT"), record("two", 20, "4M", "ACGT")]
    with pytest.raises(CallbackError, match="source"):
        StreamingPileupBuilder(broken()).pileup("chr1", 10)
    with pytest.raises(CallbackError, match="filter"):
        StreamingPileupBuilder(reads, read_filter=refuse).pileup("chr1", 10)
    with pytest.raises(CallbackError, match="filter"):
        StreamingPileupBuilder([], read_filter=refuse).accepts(reads[0])
    with pytest.raises(CallbackError, match="tap"):
        StreamingPileupBuilder(reads, tap=tap).pileup("chr1", 20)


def test_a_read_whose_cigar_and_sequence_differ_in_length_is_refused() -> None:
    reads = [record("bad", 10, "4M", "ACGTA")]
    with pytest.raises(
        ValueError, match="Read bad is invalid: CIGAR and query sequence lengths differ."
    ):
        StreamingPileupBuilder(reads).pileup("chr1", 10)


def test_builder_reads_a_cram_through_pysam(tmp_path: Path) -> None:
    contigs = {"chr1": "ACGT" * 250, "chr2": "TTGG" * 250}
    fasta = tmp_path / "reference.fa"
    fasta.write_text("".join(f">{name}\n{bases}\n" for name, bases in contigs.items()))
    pysam.faidx(str(fasta))
    cram = tmp_path / "reads.cram"
    with AlignmentFile(str(cram), "wc", header=HEADER, reference_filename=str(fasta)) as sink:
        sink.write(record("r", 10, "4M", "GTAC"))
    with (
        AlignmentFile(str(cram), reference_filename=str(fasta), threads=2) as reads,
        StreamingPileupBuilder(reads) as builder,
    ):
        assert [pileup.bases for pileup in builder.columns("chr1", 10, 14)] == [
            ["G"],
            ["T"],
            ["A"],
            ["C"],
        ]


def test_dropping_a_builder_without_closing_it_hands_no_more_reads_to_the_tap() -> None:
    reads = [record(f"r{start}", start, "4M", "ACGT") for start in (10, 20, 30)]
    tapped: list[AlignedSegment] = []
    builder = StreamingPileupBuilder(reads, tap=tapped.append)
    builder.pileup("chr1", 20)
    del builder
    assert tapped == reads[:1]


@pytest.mark.parametrize(
    "cigar,bases,placed",
    [
        ("4M", "ACGT", True),
        ("1D3M", "ACG", True),
        ("2S2N2M", "ACGT", True),
        ("1=1X2S", "ACGT", True),
        ("4D", "*", True),
        ("4S", "ACGT", False),
        ("4I", "ACGT", False),
        ("2S2I", "ACGT", False),
        ("4H4S", "ACGT", False),
        ("1S2I1S", "ACGT", False),
    ],
)
def test_a_builder_accepts_only_reads_with_a_reference_consuming_operator(
    cigar: str, bases: str, placed: bool
) -> None:
    assert StreamingPileupBuilder([]).accepts(record("r", 10, cigar, bases)) is placed


def test_a_builder_accepts_only_mapped_reads() -> None:
    builder = StreamingPileupBuilder([])
    assert not builder.accepts(unmapped("u"))
    assert not builder.accepts(record("r", 10, "4M", "ACGT", flag=4))


def test_a_builder_is_used_and_dropped_on_other_threads() -> None:
    class Sink:
        def write(self, _read: AlignedSegment) -> None:
            pass

    reads = [record("r", 10, "4M", "ACGT")]
    sink = Sink()
    alive = weakref.ref(sink)
    builders = [StreamingPileupBuilder(reads, tap=sink.write)]
    del sink
    with ThreadPoolExecutor(1) as pool:
        assert pool.submit(lambda: builders[0].pileup("chr1", 11).bases).result() == ["C"]
        columns = builders[0].columns("chr1", 12, 14)
        swept = pool.submit(list, columns).result()
        assert [pileup.bases for pileup in swept] == [["G"], ["T"]]
        del columns
        pool.submit(builders.clear).result()
    assert alive() is None
