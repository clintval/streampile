import pytest
from pysam import AlignedSegment

from streampile import AgreementStrategy
from streampile import DisagreementStrategy
from streampile import Pileup
from streampile import PileupRead
from streampile import PileupReadType
from streampile import PileupTemplate
from streampile import StreamingPileupBuilder
from streampile._pileup import BASE
from streampile._pileup import DELETION
from streampile._pileup import INSERTION
from streampile._pileup import SKIP

from .records import HEADER
from .records import entries
from .records import record


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


def called(template: PileupTemplate) -> tuple[str | None, str, str | None, int | None]:
    """A template's name, type, base, and quality."""
    return (template.query_name, template.pileup_type.value, template.base, template.qual)


def mates(base1: str, qual1: int, base2: str, qual2: int) -> Pileup:
    """The pileup at 12 of two overlapping mates of an FR pair whose bases there are given."""
    reads = [
        record("t", 10, "4M", f"AA{base1}A", flag=99, quals=[30, 30, qual1, 30]),
        record("t", 10, "4M", f"AA{base2}A", flag=147, quals=[30, 30, qual2, 30]),
    ]
    return Pileup.from_alignments(reads, "chr1", 12)


def test_templates_group_reads_by_name_in_the_order_of_their_first_entries() -> None:
    reads = [
        record("q3", 50, "50M", "A" * 50, flag=99),
        record("q1", 100, "50M", "C" * 50, flag=99),
        record("q2", 100, "50M", "G" * 50, flag=147),
        record("q3", 100, "50M", "T" * 50, flag=147),
        record("q1", 110, "50M", "C" * 50, flag=147),
        record("q2", 110, "50M", "G" * 50, flag=99),
    ]
    pileup = StreamingPileupBuilder(reads).pileup("chr1", 125)
    templates = pileup.templates()
    assert [called(template) for template in templates] == [
        ("q1", "base", "C", 80),
        ("q2", "base", "G", 80),
        ("q3", "base", "T", 40),
    ]
    assert [len(template.reads) for template in templates] == [2, 2, 1]
    assert [read.alignment.flag for read in templates[1].reads] == [147, 99]
    assert all(read in pileup.pileups for read in templates[0].reads)


@pytest.mark.parametrize(
    "agreement,quals,qual",
    [
        (AgreementStrategy.consensus, (30, 35), 65),
        (AgreementStrategy.consensus, (60, 50), 93),
        (AgreementStrategy.max_qual, (30, 35), 35),
        (AgreementStrategy.pass_through, (30, 35), 35),
    ],
)
def test_agreeing_bases_make_the_quality_of_the_agreement_strategy(
    agreement: AgreementStrategy, quals: tuple[int, int], qual: int
) -> None:
    pileup = mates("C", quals[0], "C", quals[1])
    assert called(pileup.templates(agreement=agreement)[0]) == ("t", "base", "C", qual)


@pytest.mark.parametrize(
    "disagreement,quals,base,qual",
    [
        (DisagreementStrategy.consensus, (30, 20), "A", 10),
        (DisagreementStrategy.consensus, (20, 30), "C", 10),
        (DisagreementStrategy.consensus, (21, 20), "A", 2),
        (DisagreementStrategy.consensus, (30, 30), "N", 2),
        (DisagreementStrategy.mask_both, (30, 20), "N", 2),
        (DisagreementStrategy.mask_lower_qual, (30, 20), "A", 30),
        (DisagreementStrategy.mask_lower_qual, (20, 30), "C", 30),
        (DisagreementStrategy.mask_lower_qual, (30, 30), "N", 2),
    ],
)
def test_disagreeing_bases_make_the_base_and_quality_of_the_disagreement_strategy(
    disagreement: DisagreementStrategy, quals: tuple[int, int], base: str, qual: int
) -> None:
    pileup = mates("A", quals[0], "C", quals[1])
    template = pileup.templates(disagreement=disagreement)[0]
    assert (template.base, template.qual, template.is_no_call) == (base, qual, base == "N")


def test_strategies_are_named_by_their_values() -> None:
    pileup = mates("A", 30, "C", 20)
    named = pileup.templates(
        agreement="max_qual",  # type: ignore[arg-type]  # pyright: ignore[reportArgumentType]  # ty: ignore[invalid-argument-type]
        disagreement="mask_lower_qual",  # type: ignore[arg-type]  # pyright: ignore[reportArgumentType]  # ty: ignore[invalid-argument-type]
    )
    assert (named[0].base, named[0].qual) == ("A", 30)
    with pytest.raises(ValueError, match="'bogus' is not a valid AgreementStrategy"):
        pileup.templates(agreement="bogus")  # type: ignore[arg-type]  # pyright: ignore[reportArgumentType]  # ty: ignore[invalid-argument-type]
    with pytest.raises(ValueError, match="'bogus' is not a valid DisagreementStrategy"):
        pileup.templates(disagreement="bogus")  # type: ignore[arg-type]  # pyright: ignore[reportArgumentType]  # ty: ignore[invalid-argument-type]


def test_each_unnamed_read_is_a_template_named_star_of_its_own() -> None:
    reads = [
        AlignedSegment.fromstring(f"*\t0\tchr1\t11\t60\t4M\t*\t0\t0\t{bases}\tIIII", HEADER)
        for bases in ("AAAA", "AAAA", "CCCC")
    ]
    pileup = StreamingPileupBuilder(reads).pileup("chr1", 11)
    assert [called(template) for template in pileup.templates()] == [
        ("*", "base", "A", 40),
        ("*", "base", "A", 40),
        ("*", "base", "C", 40),
    ]


def test_a_no_call_leaves_the_other_reads_base_at_its_own_quality() -> None:
    assert called(mates("N", 40, "A", 20).templates()[0]) == ("t", "base", "A", 20)
    assert called(mates("N", 10, "N", 30).templates()[0]) == ("t", "base", "N", 30)


def test_a_read_with_a_deletion_or_a_skip_holds_no_base() -> None:
    base = record("t", 10, "4M", "ACGT", flag=99, quals=[30, 30, 15, 30])
    deletion = record("t", 10, "2M1D2M", "ACTT", flag=147, quals=[30, 30, 25, 30])
    skip = record("t", 10, "2M1N2M", "ACTT", flag=147)
    no_call = record("t", 10, "4M", "ACNT", flag=99)
    for reads, expected in (
        ([base, deletion], ("t", "base", "G", 15)),
        ([no_call, deletion], ("t", "base", "N", 40)),
        ([base, skip], ("t", "base", "G", 15)),
        ([deletion], ("t", "deletion", None, 25)),
        ([skip], ("t", "skip", None, None)),
    ):
        template = Pileup.from_alignments(reads, "chr1", 12).templates()[0]
        assert called(template) == expected
        assert (template.is_del, template.is_refskip) == (
            expected[1] == "deletion",
            expected[1] == "skip",
        )


def test_insertion_entries_are_no_part_of_a_template() -> None:
    reads = [
        record("t", 10, "3M2I3M", "ACGTTACG", flag=99),
        record("t", 10, "6M", "ACGACG", flag=147),
        record("opens", 13, "1I3M", "TACG"),
    ]
    pileup = Pileup.from_alignments(reads, "chr1", 12)
    templates = pileup.templates()
    assert (len(pileup.pileups), [called(template) for template in templates]) == (
        4,
        [("t", "base", "G", 80)],
    )
    assert [read.pileup_type for read in templates[0].reads] == [BASE, BASE]


def test_a_read_under_the_floor_does_not_vote() -> None:
    pileup = mates("C", 10, "C", 10)
    assert (pileup.min_base_quality, pileup.filtered_depth) == (13, 0)
    template = pileup.templates()[0]
    assert (len(template.reads), called(template)) == (2, ("t", "base", None, None))
    unfloored = Pileup.from_alignments([read.alignment for read in pileup.pileups], "chr1", 12, 0)
    assert called(unfloored.templates()[0]) == ("t", "base", "C", 20)


def test_a_mate_under_the_floor_does_not_mask_the_other_mates_base() -> None:
    pileup = mates("A", 35, "C", 5)
    masked = pileup.templates(disagreement=DisagreementStrategy.mask_both)[0]
    assert (masked.base, masked.qual) == ("A", 35)
    assert (pileup.templates()[0].base, pileup.templates()[0].qual) == ("A", 35)


def test_a_templates_strand_and_distances_are_its_first_reads() -> None:
    reads = [
        record("t", 100, "10M", "A" * 10, flag=99),
        record("t", 105, "10M", "C" * 10, flag=147),
    ]
    for read, mate in ((reads[0], reads[1]), (reads[1], reads[0])):
        read.next_reference_id = 0
        read.next_reference_start = mate.reference_start
        read.set_tag("MC", "10M")  # pyright: ignore[reportUnknownMemberType]
    distances: list[tuple[bool, int | None, int | None]] = []
    for pos in (100, 105, 110, 114):
        template = Pileup.from_alignments(reads, "chr1", pos).templates()[0]
        distances.append((
            template.is_reverse,
            template.five_prime_distance,
            template.template_end_distance,
        ))
    assert distances == [(False, 0, 14), (False, 5, 9), (False, 10, 4), (False, 14, 0)]
    reads[1].set_tag("MC", None)  # pyright: ignore[reportUnknownMemberType]
    template = Pileup.from_alignments(reads, "chr1", 110).templates()[0]
    with pytest.raises(ValueError, match="no MC tag"):
        _ = template.five_prime_distance
    assert template.template_end_distance == 4


def test_a_template_has_both_distances_where_it_holds_a_deletion_or_a_skip() -> None:
    def mated(cigar1: str, start2: int, cigar2: str) -> list[AlignedSegment]:
        reads = [
            record("t", 999, cigar1, "A" * 100, flag=99),
            record("t", start2, cigar2, "C" * 100, flag=147),
        ]
        for read, mate in ((reads[0], reads[1]), (reads[1], reads[0])):
            read.next_reference_id = 0
            read.next_reference_start = mate.reference_start
            read.set_tag("MC", mate.cigarstring)  # pyright: ignore[reportUnknownMemberType]
        return reads

    def distances(reads: list[AlignedSegment], pos: int) -> tuple[int | None, int | None]:
        template = Pileup.from_alignments(reads, "chr1", pos).templates()[0]
        assert template.pileup_type in (DELETION, SKIP)
        return template.five_prime_distance, template.template_end_distance

    deleted = mated("50M2D50M", 1199, "100M")
    assert [distances(deleted, pos) for pos in (1049, 1050)] == [(50, 248)] * 2
    skipped = mated("50M100N50M", 1299, "100M")
    assert [distances(skipped, pos) for pos in (1049, 1148)] == [(50, 250)] * 2
    deleted_second = mated("100M", 1099, "50M2D50M")
    assert distances(deleted_second, 1149) == (150, 50)
    (entry,) = Pileup.from_alignments(deleted_second, "chr1", 1149).pileups
    assert (entry.five_prime_distance, entry.template_end_distance) == (50, 150)


def test_a_read_is_of_an_fr_pair_when_its_forward_5_prime_end_is_at_or_before_its_reverse() -> None:
    def mated(start1: int, flag1: int, start2: int, flag2: int) -> list[AlignedSegment]:
        reads = [
            record("t", start1, "10M", "A" * 10, flag=flag1),
            record("t", start2, "10M", "C" * 10, flag=flag2),
        ]
        for read, mate in ((reads[0], reads[1]), (reads[1], reads[0])):
            read.next_reference_id = 0
            read.next_reference_start = mate.reference_start
            read.set_tag("MC", "10M")  # pyright: ignore[reportUnknownMemberType]
        return reads

    tie = mated(109, 99, 100, 147)
    pileup = Pileup.from_alignments(tie, "chr1", 109)
    assert [(read.is_fr_pair, read.template_end_distance) for read in pileup.pileups] == [
        (True, 0),
        (True, 0),
    ]
    outward = mated(100, 83, 200, 163)
    for pos in (105, 205):
        (read,) = Pileup.from_alignments(outward, "chr1", pos).pileups
        assert (read.is_fr_pair, read.template_end_distance) == (False, None)
    tie[0].set_tag("MC", None)  # pyright: ignore[reportUnknownMemberType]
    (read, _) = Pileup.from_alignments(tie, "chr1", 109).pileups
    with pytest.raises(ValueError, match="no MC tag"):
        _ = read.is_fr_pair


def test_a_template_reads_its_strand_from_its_second_read_without_its_first() -> None:
    read = record("t", 10, "4M", "ACGT", flag=163)
    template = Pileup.from_alignments([read], "chr1", 12).templates()[0]
    assert (template.is_reverse, template.five_prime_distance) == (True, None)
    assert repr(template).startswith("PileupTemplate(query_name='t', pileup_type=")


def test_from_alignments_drops_reads_on_other_contigs() -> None:
    reads = [record("a", 10, "4M", "ACGT"), record("b", 10, "4M", "TTTT", contig="chr2")]
    for contig, name, base in (("chr1", "a", "C"), ("chr2", "b", "T")):
        pileup = Pileup.from_alignments(reads, contig, 11)
        assert [entry.alignment.query_name for entry in pileup.pileups] == [name]
        assert pileup.bases == [base]


def test_a_pileup_is_frozen_and_compared_and_hashed_by_its_fields() -> None:
    read = record("r", 10, "4M", "ACGT")
    pileup = Pileup("chr1", 11, (PileupRead(read, 1, 1, BASE),))
    same = Pileup("chr1", 11, [PileupRead(read, 1, 1, BASE)], min_base_quality=13)
    assert (pileup, hash(pileup)) == (same, hash(same))
    assert pileup != Pileup("chr1", 11, (PileupRead(read, 1, 1, BASE),), min_base_quality=14)
    assert repr(pileup).startswith("Pileup(reference_name='chr1', reference_pos=11, pileups=(")
    with pytest.raises(AttributeError):
        # Each checker fails on an unused ignore, so all three must reject this assignment.
        pileup.reference_pos = 12  # type: ignore[misc]  # pyright: ignore[reportAttributeAccessIssue]  # ty: ignore[invalid-assignment]


def test_a_pileup_read_behaves_as_the_tuple_of_its_fields() -> None:
    read = record("r", 10, "2M2I2M", "ACGTAC")
    entry = PileupRead(read, None, None, INSERTION, insertion_offset=2, insertion_length=2)
    fields = (read, None, None, INSERTION, 2, 2)
    assert (entry == fields, tuple(entry), len(entry), entry[3], entry[-1]) == (
        True,
        fields,
        6,
        INSERTION,
        2,
    )
    assert hash(entry) == hash(fields)
    assert entry._fields == (
        "alignment",
        "query_position",
        "query_position_or_next",
        "pileup_type",
        "insertion_offset",
        "insertion_length",
    )
    assert entry._asdict() == dict(zip(entry._fields, fields, strict=True))
    moved = entry._replace(insertion_offset=3, insertion_length=1)
    assert (moved.inserted_bases, moved.alignment) == ("T", read)
    assert entry < moved and moved != entry
    assert repr(entry).startswith("PileupRead(alignment=")
    with pytest.raises(ValueError, match="Got unexpected field names"):
        entry._replace(base="A")  # type: ignore[call-arg]  # pyright: ignore[reportCallIssue]  # ty: ignore[unknown-argument]
    with pytest.raises(ValueError, match="'bogus' is not a valid PileupReadType"):
        PileupRead(read, 0, 0, "bogus")  # type: ignore[arg-type]  # pyright: ignore[reportArgumentType]  # ty: ignore[invalid-argument-type]
