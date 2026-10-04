from array import array

import pytest
from pysam import AlignedSegment

from streampile import Pileup
from streampile import PileupRead
from streampile import PileupReadType
from streampile import StreamingPileupBuilder

from .records import HEADER
from .records import entries
from .records import record

BASE = PileupReadType.base
DELETION = PileupReadType.deletion
INSERTION = PileupReadType.insertion
SKIP = PileupReadType.skip


def columns(reads: list[AlignedSegment], positions: range, **options: int) -> list[Pileup]:
    """The pileup of the reads at each position, built from the reads alone."""
    return [Pileup.from_alignments(reads, "chr1", pos, **options) for pos in positions]


def test_a_base_entry_holds_its_base_and_quality() -> None:
    read = record("r", 10, "4M", "ACGT", quals=[10, 20, 30, 40])
    entry = PileupRead(read, 2, 2, BASE)
    assert (entry.base, entry.qual) == ("G", 30)


def test_an_insertion_entry_holds_no_base_or_quality() -> None:
    read = record("r", 10, "2M2I2M", "ACGTAC", quals=[10, 20, 30, 40, 50, 60])
    entry = PileupRead(read, None, None, INSERTION, insertion_offset=2, insertion_length=2)
    assert (entry.base, entry.qual) == (None, None)
    assert (entry.inserted_bases, entry.inserted_qualities) == ("GT", [30, 40])


def test_a_deletion_entry_holds_the_quality_of_the_next_base() -> None:
    read = record("r", 10, "1M1D3M", "ACGT", quals=[10, 20, 30, 40])
    entry = PileupRead(read, None, 1, DELETION)
    assert (entry.base, entry.qual) == (None, 20)
    assert (entry.inserted_bases, entry.inserted_qualities) == (None, None)


def test_a_base_entry_of_an_n_is_a_no_call() -> None:
    read = record("r", 10, "4M", "AnGT")
    entries = [PileupRead(read, offset, offset, BASE) for offset in range(4)]
    assert [entry.is_no_call for entry in entries] == [False, True, False, False]
    assert not PileupRead(read, None, None, SKIP).is_no_call


@pytest.mark.parametrize(
    "pileup_type,is_del,is_ins,is_refskip",
    [
        (BASE, False, False, False),
        (DELETION, True, False, False),
        (INSERTION, False, True, False),
        (SKIP, False, False, True),
    ],
)
def test_an_entry_is_one_kind(
    pileup_type: PileupReadType, is_del: bool, is_ins: bool, is_refskip: bool
) -> None:
    entry = PileupRead(record("r", 10, "4M", "ACGT"), None, None, pileup_type)
    assert (entry.is_del, entry.is_ins, entry.is_refskip) == (is_del, is_ins, is_refskip)


def test_depths_count_bases_deletions_and_skips_but_not_insertions() -> None:
    base = PileupRead(record("base", 5, "10M", "A" * 10), 5, 5, BASE)
    deletion = PileupRead(record("deletion", 5, "3M3D4M", "AAATTTT"), None, 3, DELETION)
    skip = PileupRead(record("skip", 5, "3M3N4M", "AAATTTT"), None, None, SKIP)
    insertion = PileupRead(record("insertion", 11, "1I6M", "AAATTTT"), None, None, INSERTION, 0, 1)
    depths = {
        name: (pileup.unfiltered_depth, pileup.filtered_depth)
        for name, pileup in {
            "empty": Pileup("chr1", 10, ()),
            "base": Pileup("chr1", 10, (base,)),
            "deletion": Pileup("chr1", 10, (deletion,)),
            "skip": Pileup("chr1", 10, (skip,)),
            "insertion": Pileup("chr1", 10, (insertion,)),
            "mixed": Pileup("chr1", 10, (base, deletion, skip, insertion)),
        }.items()
    }
    assert depths == {
        "empty": (0, 0),
        "base": (1, 1),
        "deletion": (1, 1),
        "skip": (1, 0),
        "insertion": (0, 0),
        "mixed": (3, 2),
    }


def test_qualities_of_one_read() -> None:
    reads = [record("r", 1, "4M", "ACGT", quals=[20, 21, 22, 23])]
    assert [pileup.qualities for pileup in columns(reads, range(1, 5))] == [
        [20],
        [21],
        [22],
        [23],
    ]
    floored = columns(reads, range(1, 5), min_base_quality=22)
    assert [pileup.qualities for pileup in floored] == [[], [], [22], [23]]


def test_qualities_leave_out_deletions_and_insertions() -> None:
    deleted = columns([record("r", 1, "2M1D1M", "ACG", quals=[20, 21, 22])], range(1, 5))
    assert [pileup.qualities for pileup in deleted] == [[20], [21], [], [22]]
    assert deleted[2].filtered_depth == 1
    inserted = columns([record("r", 1, "1M3I1M", "AGGGT", quals=[31, 32, 33, 34, 35])], range(1, 4))
    assert [pileup.qualities for pileup in inserted] == [[31], [35], []]
    assert entries(inserted[0])[1] == ("r", "insertion", None, None, "GGG")


def test_qualities_of_overlapping_reads() -> None:
    reads = [
        record("one", 1, "4M", "ACGT", quals=[20, 21, 22, 23]),
        record("two", 3, "4M", "TGCA", quals=[30, 31, 32, 33]),
    ]
    assert [pileup.qualities for pileup in columns(reads, range(3, 5))] == [
        [22, 30],
        [23, 31],
    ]


def test_bases_of_one_read() -> None:
    reads = [record("r", 1, "4M", "ACGT", quals=[12, 13, 25, 30])]
    assert [pileup.bases for pileup in columns(reads, range(0, 6))] == [
        [],
        [],
        ["C"],
        ["G"],
        ["T"],
        [],
    ]
    unfloored = columns(reads, range(1, 2), min_base_quality=0)
    assert unfloored[0].bases == ["A"]


def test_bases_leave_out_deletions_and_insertions() -> None:
    deleted = columns([record("r", 1, "2M1D1M", "ACG")], range(1, 5))
    assert [pileup.bases for pileup in deleted] == [["A"], ["C"], [], ["G"]]
    inserted = columns([record("r", 1, "1M3I1M", "AGGGT")], range(1, 4))
    assert [pileup.bases for pileup in inserted] == [["A"], ["T"], []]


def test_bases_of_overlapping_reads() -> None:
    reads = [record("one", 1, "4M", "ACGT"), record("two", 3, "4M", "TGCA")]
    assert [pileup.bases for pileup in columns(reads, range(3, 5))] == [
        ["G", "T"],
        ["T", "G"],
    ]


def test_views_of_reads_with_no_stored_bases_or_no_cigar_are_empty() -> None:
    no_bases = columns([record("r", 0, "4M", "*")], range(0, 1), min_base_quality=0)[0]
    assert (no_bases.unfiltered_depth, no_bases.filtered_depth) == (1, 0)
    assert (no_bases.bases, no_bases.qualities) == ([], [])
    no_cigar = AlignedSegment(HEADER)
    no_cigar.reference_name = "chr1"
    no_cigar.reference_start = 0
    no_cigar.query_sequence = None
    unplaced = Pileup.from_alignments([no_cigar], "chr1", 0)
    assert (unplaced.pileups, unplaced.bases, unplaced.qualities) == ((), [], [])


def test_without_overlaps_keeps_the_first_read_of_each_template() -> None:
    reads = [
        record("q3", 50, "50M", "A" * 50, flag=99),
        record("q1", 100, "50M", "C" * 50, flag=99),
        record("q2", 100, "50M", "G" * 50, flag=147),
        record("q3", 100, "50M", "T" * 50, flag=147),
        record("q1", 110, "50M", "C" * 50, flag=147),
        record("q2", 110, "50M", "G" * 50, flag=99),
    ]
    pileup = StreamingPileupBuilder(reads).pileup("chr1", 125)
    kept = pileup.without_overlaps()
    assert (pileup.unfiltered_depth, kept.unfiltered_depth) == (5, 3)
    assert [(entry.alignment.query_name, entry.alignment.flag) for entry in kept.pileups] == [
        ("q1", 99),
        ("q2", 147),
        ("q3", 147),
    ]
    assert kept.bases == ["C", "G", "T"]


def test_without_overlaps_keeps_every_entry_of_the_kept_read() -> None:
    reads = [
        record("pair", 10, "3M2I3M", "ACGTTACG", flag=99),
        record("pair", 10, "6M", "ACGACG", flag=147),
        record("other", 12, "4M", "GACG"),
    ]
    pileup = Pileup.from_alignments(reads, "chr1", 12, min_base_quality=30)
    kept = pileup.without_overlaps()
    assert entries(kept) == [
        ("pair", "base", 2, 2, None),
        ("pair", "insertion", None, None, "TT"),
        ("other", "base", 0, 0, None),
    ]
    assert [entry.alignment.flag for entry in kept.pileups] == [99, 99, 0]
    assert (kept.reference_name, kept.reference_pos, kept.min_base_quality) == ("chr1", 12, 30)
    assert len(pileup.pileups) == 4


def test_without_overlaps_keeps_a_mate_whose_base_another_mate_skips() -> None:
    reads = [
        record("t", 100, "20M300N20M", "A" * 40, flag=99),
        record("t", 330, "40M", "G" * 10 + "C" + "G" * 29, flag=147),
    ]
    pileup = Pileup.from_alignments(reads, "chr1", 340)
    kept = pileup.without_overlaps()
    assert [entry.pileup_type for entry in pileup.pileups] == [SKIP, BASE]
    assert (kept.filtered_depth, kept.bases) == (1, ["C"])
    assert [entry.alignment.flag for entry in kept.pileups] == [147]


def test_without_overlaps_keeps_a_mate_at_the_floor_over_one_under_it() -> None:
    reads = [
        record("t", 100, "10M", "A" * 10, flag=99, quals=[2] * 10),
        record("t", 105, "10M", "GGC" + "G" * 7, flag=147),
    ]
    kept = Pileup.from_alignments(reads, "chr1", 107, min_base_quality=13).without_overlaps()
    assert (kept.filtered_depth, kept.bases) == (1, ["C"])
    reads[1].query_qualities = array("B", [2] * 10)
    kept = Pileup.from_alignments(reads, "chr1", 107, min_base_quality=13).without_overlaps()
    assert [entry.alignment.flag for entry in kept.pileups] == [99]
    assert kept.filtered_depth == 0


def test_from_alignments_drops_reads_on_other_contigs() -> None:
    reads = [record("a", 10, "4M", "ACGT"), record("b", 10, "4M", "TTTT", contig="chr2")]
    for contig, name, base in (("chr1", "a", "C"), ("chr2", "b", "T")):
        pileup = Pileup.from_alignments(reads, contig, 11)
        assert [entry.alignment.query_name for entry in pileup.pileups] == [name]
        assert pileup.bases == [base]
