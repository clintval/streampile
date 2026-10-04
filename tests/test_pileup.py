from streampile import Pileup
from streampile import StreamingPileupBuilder

from .records import entries
from .records import record


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
    assert kept.get_query_sequences == ["C", "G", "T"]


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
