import random
from collections import Counter
from pathlib import Path

import pytest
from pysam import AlignedSegment
from pysam import AlignmentFile
from pysam import AlignmentHeader
from pysam import FastaFile

from streampile import StreamingPileupBuilder
from streampile import TabulatedBase
from streampile import TabulationReader
from streampile import TabulationWriter
from streampile import Tabulator
from streampile import normalize
from streampile import tabulate

from .records import DATA
from .records import header_of
from .records import record
from .records import territory
from .records import unmapped
from .records import write_bam
from .records import write_fasta

CHR1 = "ACGTAGGCTAACGTTAGCCATGCAAAAAGTCCATGACGTCGATCGGATCCTAGGCTAGCT"
HEADER = AlignmentHeader.from_text(
    "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:60\n@SQ\tSN:chr2\tLN:40\n"
)


@pytest.fixture
def reference() -> FastaFile:
    return FastaFile(str(DATA / "reference.fa"))


def alleles(base: TabulatedBase) -> dict[str, int]:
    return {
        f"{ref}>{alt}": reads
        for ref, alt, reads in zip(base.alt_refs, base.alts, base.alt_reads, strict=True)
    }


def sites_of(
    spans: list[tuple[str, int, int]], **options: int
) -> dict[tuple[str, int], TabulatedBase]:
    with (
        AlignmentFile(str(DATA / "reads.bam")) as reads,
        FastaFile(str(DATA / "reference.fa")) as fasta,
    ):
        return {
            (site.contig, site.pos): site
            for site in tabulate(reads, fasta, territory(*spans), **options)
        }


def tabulated(
    tmp_path: Path, chr1: str, reads: list[AlignedSegment], **options: int
) -> list[TabulatedBase]:
    """Tabulate reads, made with `header_of({"chr1": chr1})`, over the whole of `chr1`."""
    fasta = write_fasta(tmp_path / "reference.fa", {"chr1": chr1})
    path = write_bam(tmp_path / "reads.bam", reads, header=header_of({"chr1": chr1}))
    with AlignmentFile(str(path)) as alignments, FastaFile(str(fasta)) as reference:
        return list(tabulate(alignments, reference, territory(("chr1", 0, len(chr1))), **options))


def spanning(sites: list[TabulatedBase]) -> list[int]:
    """The reads at each base with an allele anchored at an earlier base that spans it."""
    counts = Counter(
        site.pos + offset
        for site in sites
        for ref, reads in zip(site.alt_refs, site.alt_reads, strict=True)
        for offset in range(1, len(ref))
        for _ in range(reads)
    )
    return [counts[site.pos] for site in sites]


def alleles_of(
    tabulator: Tabulator, start: int, cigar: str, bases: str, quals: list[int] | None = None
) -> list[tuple[int, str]]:
    counted, _ = tabulator.alleles(record("r", start, cigar, bases, quals=quals, header=HEADER))
    return [(allele.pos, allele.key) for allele in counted]


def test_the_reference_matches_the_fixture(reference: FastaFile) -> None:
    assert reference.fetch("chr1") == CHR1


@pytest.mark.parametrize(
    "pos,ref,alt,expected",
    [
        (10, "A", "T", (10, "A", "T")),
        (10, "AC", "GG", (10, "AC", "GG")),
        (26, "AA", "A", (22, "CA", "C")),
        (33, "T", "TTT", (32, "A", "ATT")),
        (10, "ACG", "AG", (10, "AC", "A")),
        (10, "ACGT", "ACCT", (12, "G", "C")),
        (25, "AAA", "AA", None),
    ],
)
def test_normalize_trims_and_left_aligns(
    pos: int, ref: str, alt: str, expected: tuple[int, str, str] | None
) -> None:
    floor = 23 if expected is None else 0
    assert normalize(pos, ref, alt, lambda start, end: CHR1[start:end], floor) == expected


def test_normalize_returns_at_once_for_an_allele_that_changes_nothing() -> None:
    asked: list[tuple[int, int]] = []

    def reference(start: int, end: int) -> str:
        asked.append((start, end))
        return "A" * (end - start)

    assert normalize(1_000_000, "A", "A", reference) is None
    assert normalize(1_000_000, "AC", "AC", reference) is None
    assert asked == []


def test_alleles_of_one_read(reference: FastaFile) -> None:
    tabulator = Tabulator(reference, min_base_quality=30)
    assert alleles_of(tabulator, 0, "20M", CHR1[0:20]) == []
    assert alleles_of(tabulator, 8, "4M", "TTGC") == [(9, "AA>TG")]
    assert alleles_of(tabulator, 8, "4M", "TACC") == [(10, "A>C")]
    assert alleles_of(tabulator, 20, "7M1D3M", CHR1[20:27] + CHR1[28:31]) == [(22, "CA>C")]
    assert alleles_of(tabulator, 26, "8M2I2M", CHR1[26:34] + "TT" + CHR1[34:36]) == [(32, "A>ATT")]
    assert alleles_of(tabulator, 10, "2M1D2M", "AATT") == [(11, "CG>A")]
    assert alleles_of(tabulator, 10, "2M2I2M", "ACGGGT") == [(11, "C>CGG")]


def test_reads_are_not_counted_for_alleles_they_cannot_place(reference: FastaFile) -> None:
    tabulator = Tabulator(reference, min_base_quality=30)
    assert alleles_of(tabulator, 10, "2I4M", "TTACGT") == []
    assert alleles_of(tabulator, 10, "2S1D4M", "GGCGTT") == []
    assert alleles_of(tabulator, 10, "4M2D", "ACGT") == []
    assert alleles_of(tabulator, 24, "2M1D4M", "AAAGTC") == []
    assert alleles_of(tabulator, 10, "4M", "TNGT") == []
    assert alleles_of(tabulator, 10, "4M", "TCGT", quals=[20, 40, 40, 40]) == []
    assert alleles_of(tabulator, 10, "4M2N4M", "TCGT" + CHR1[16:20]) == [(10, "A>T")]
    counted, dropped = tabulator.alleles(
        record("r", 10, "4M", "TCGT", quals=[20, 40, 40, 40], header=HEADER)
    )
    assert counted == [] and dropped == [(10, 11)]


def test_alleles_need_a_mapped_read(reference: FastaFile) -> None:
    with pytest.raises(ValueError, match="Read u is not mapped."):
        Tabulator(reference).alleles(unmapped("u", header=HEADER))


def test_tabulate_the_fixture_by_hand() -> None:
    sites = sites_of([("chr1", 0, 60)], min_base_quality=30, min_mapping_quality=20)
    assert len(sites) == 60
    assert sites["chr1", 1] == TabulatedBase(contig="chr1", pos=1, ref="A", depth=1, ref_reads=1)
    assert sites["chr1", 11] == TabulatedBase(
        contig="chr1",
        pos=11,
        ref="A",
        depth=4,
        ref_reads=2,
        alt_refs=("A", "AC"),
        alts=("T", "GG"),
        alt_reads=(1, 1),
    )
    assert (sites["chr1", 12].depth, sites["chr1", 12].ref_reads) == (4, 3)
    assert (sites["chr1", 13].depth, sites["chr1", 13].ref_reads) == (3, 3)
    assert alleles(sites["chr1", 23]) == {"CA>C": 1}
    assert (sites["chr1", 24].depth, sites["chr1", 24].ref_reads) == (4, 3)
    assert alleles(sites["chr1", 33]) == {"A>ATT": 1}
    assert [sites["chr1", pos].depth for pos in range(34, 41)] == [4, 4, 3, 2, 2, 1, 1]
    assert [sites["chr1", pos].depth for pos in range(45, 51)] == [1, 3, 3, 3, 3, 3]


def test_tabulate_counts_every_informative_read_once() -> None:
    sites = sites_of([("chr1", 0, 60), ("chr2", 0, 40)], min_base_quality=30)
    spanning = {("chr1", 12): 1, ("chr1", 24): 1}
    for key, site in sites.items():
        assert site.depth == site.ref_reads + sum(site.alt_reads) + spanning.get(key, 0)


def test_tabulate_leaves_out_reads_by_flag_and_mapping_quality() -> None:
    with_duplicates = sites_of([("chr1", 32, 42)], exclude_flags=0)
    without = sites_of([("chr1", 32, 42)])
    assert (
        sum(site.depth for site in with_duplicates.values())
        == sum(site.depth for site in without.values()) + 10
    )
    assert (
        sum(site.depth for site in sites_of([("chr1", 0, 60)], min_mapping_quality=61).values())
        == 0
    )


def test_tabulate_joins_spans_and_writes_zero_depth_bases_in_header_order() -> None:
    sites = sites_of([("chr2", 30, 38), ("chr2", 0, 5), ("chr2", 3, 8), ("chr1", 58, 60)])
    assert list(sites) == [("chr1", 59), ("chr1", 60)] + [("chr2", pos) for pos in range(1, 9)] + [
        ("chr2", pos) for pos in range(31, 39)
    ]
    assert [sites["chr2", pos].depth for pos in range(4, 9)] == [0, 0, 1, 1, 1]


def test_tabulate_refuses_a_contig_not_in_the_header_and_skips_empty_spans() -> None:
    with pytest.raises(ValueError, match="Contig chr3 is not in the alignment header."):
        sites_of([("chr1", 0, 10), ("chr3", 0, 10)])
    assert sites_of([("chr1", 10, 10)]) == {}


def test_tabulate_checks_the_territory_when_called(tmp_path: Path) -> None:
    short = write_fasta(tmp_path / "short.fa", {"chr1": CHR1[:50], "chr2": "A" * 40})
    only_chr1 = write_fasta(tmp_path / "chr1.fa", {"chr1": CHR1})
    refused = [
        (DATA / "reference.fa", ("chr3", 0, 10), "Contig chr3 is not in the alignment header."),
        (only_chr1, ("chr2", 0, 10), "Contig chr2 is not in the reference."),
        (short, ("chr1", 0, 10), "Contig chr1 has 50 bases in the reference but 60 in the"),
        (DATA / "reference.fa", ("chr2", 30, 45), "Span chr2:30-45 runs past the end of chr2,"),
    ]
    with AlignmentFile(str(DATA / "reads.bam")) as reads:
        for path, span, message in refused:
            with FastaFile(str(path)) as fasta, pytest.raises(ValueError, match=message):
                _ = tabulate(reads, fasta, territory(("chr2", 0, 5), span))


def test_tabulate_counts_reads_spanning_several_spans(tmp_path: Path, reference: FastaFile) -> None:
    reads = [record("long", 0, "30M", CHR1[0:10] + "T" + CHR1[11:30], header=HEADER)]
    path = write_bam(tmp_path / "reads.bam", reads, header=HEADER)
    with AlignmentFile(str(path)) as alignments:
        sites = list(
            tabulate(
                alignments, reference, territory(("chr1", 0, 5), ("chr1", 8, 12), ("chr1", 25, 28))
            )
        )
    assert [site.depth for site in sites] == [1] * 12
    assert [site.pos for site in sites if site.alts] == [11]


def test_tabulate_counts_bases_of_any_quality_but_not_qc_fail_reads_by_default(
    tmp_path: Path, reference: FastaFile
) -> None:
    reads = [
        record("q0", 0, "4M", CHR1[0:4], quals=[0] * 4, header=HEADER),
        record("qcfail", 0, "4M", CHR1[0:4], flag=512, header=HEADER),
    ]
    path = write_bam(tmp_path / "reads.bam", reads, header=HEADER)
    with AlignmentFile(str(path)) as alignments:
        sites = list(tabulate(alignments, reference, territory(("chr1", 0, 4))))
    assert [(site.depth, site.ref_reads) for site in sites] == [(1, 1)] * 4


def test_tabulate_reads_base_qualities_of_95_and_more(tmp_path: Path, reference: FastaFile) -> None:
    quals = [40, 40, 96, 40, 120, 40, 10, 40]
    reads = [record("r", 0, "8M", CHR1[0:2] + "T" + CHR1[3:8], quals=quals, header=HEADER)]
    path = write_bam(tmp_path / "reads.bam", reads, header=HEADER)
    with AlignmentFile(str(path)) as alignments:
        sites = list(
            tabulate(alignments, reference, territory(("chr1", 0, 8)), min_base_quality=30)
        )
        assert StreamingPileupBuilder(alignments.fetch("chr1")).pileup("chr1", 2).qualities == [96]
    assert [site.depth for site in sites] == [1, 1, 1, 1, 1, 1, 0, 1]
    assert [(site.pos, alleles(site)) for site in sites if site.alts] == [(3, {"G>T": 1})]


def test_tabulate_counts_reads_with_no_stored_qualities_at_every_floor(
    tmp_path: Path, reference: FastaFile
) -> None:
    reads = [
        record("noquals", 0, "10M", CHR1[0:5] + "T" + CHR1[6:10], quals="*", header=HEADER),
        record("nobases", 0, "10M", "*", header=HEADER),
    ]
    path = write_bam(tmp_path / "reads.bam", reads, header=HEADER)
    with AlignmentFile(str(path)) as alignments:
        sites = list(
            tabulate(alignments, reference, territory(("chr1", 0, 10)), min_base_quality=60)
        )
    assert [site.depth for site in sites] == [1] * 10
    assert [(site.pos, alleles(site)) for site in sites if site.alts] == [(6, {"G>T": 1})]


def test_tabulate_agrees_with_counting_bases_in_pileups(
    tmp_path: Path, reference: FastaFile
) -> None:
    reads = [
        record("a", 0, "25M", CHR1[0:7] + "T" + CHR1[8:25], header=HEADER),
        record(
            "b",
            3,
            "25M",
            CHR1[3:7] + "T" + CHR1[8:20] + "G" + CHR1[21:28],
            quals=[40] * 17 + [10] + [40] * 7,
            header=HEADER,
        ),
        record("c", 5, "2S20M", "GG" + CHR1[5:15] + "C" + CHR1[16:25], header=HEADER),
        record("d", 10, "25M", CHR1[10:35], flag=1024, header=HEADER),
    ]
    path = write_bam(tmp_path / "reads.bam", reads, header=HEADER)
    with AlignmentFile(str(path)) as alignments:
        sites = list(
            tabulate(alignments, reference, territory(("chr1", 0, 40)), min_base_quality=30)
        )
    with (
        AlignmentFile(str(path)) as alignments,
        StreamingPileupBuilder(alignments, min_base_quality=30) as builder,
    ):
        columns = list(builder.columns("chr1", 0, 40))
    for site, column in zip(sites, columns, strict=True):
        bases = Counter(column.bases)
        assert site.depth == bases.total()
        assert site.ref_reads == bases[site.ref]
        assert alleles(site) == {
            f"{site.ref}>{base}": count for base, count in bases.items() if base != site.ref
        }


def test_bases_round_trip_through_a_table(tmp_path: Path) -> None:
    bases = list(sites_of([("chr1", 8, 40)], min_base_quality=30).values())
    with TabulationWriter.from_path(tmp_path / "bases.tsv.gz") as writer:
        writer.write_header()
        for base in bases:
            writer.write(base)
    assert list(TabulationReader.from_path(tmp_path / "bases.tsv.gz")) == bases
    assert any(base.alts for base in bases) and any(not base.alts for base in bases)


def test_left_alignment_stops_at_the_end_of_a_reference_skip(tmp_path: Path) -> None:
    chr1 = "ACGT" + "A" * 10 + "CGTACGTACG"
    reads = [record("r", 0, "4M4N2M1D3M", "ACGTAAAAA", header=header_of({"chr1": chr1}))]
    sites = tabulated(tmp_path, chr1, reads)
    assert [site.depth for site in sites[:14]] == [1] * 4 + [0] * 7 + [1] * 3
    assert not any(site.alts for site in sites)


def test_left_alignment_stops_at_the_end_of_the_previous_difference(tmp_path: Path) -> None:
    chr1 = "TTC" + "AAAA" + "GCATGCATGC"
    reads = [record("r", 0, "6M1D4M", "TTTAAA" + "GCAT", header=header_of({"chr1": chr1}))]
    sites = tabulated(tmp_path, chr1, reads)
    assert [site.depth for site in sites[:11]] == [1] * 3 + [0] * 4 + [1] * 4
    assert [(site.pos, alleles(site)) for site in sites if site.alts] == [(3, {"C>T": 1})]
    assert [site.ref_reads for site in sites[:3]] == [1, 1, 0]


def test_left_alignment_stops_at_a_read_n(tmp_path: Path) -> None:
    chr1 = "GGCAAAAGTCATGCA"
    reads = [record("r", 0, "6M1D5M", "GGCNAA" + "GTCAT", header=header_of({"chr1": chr1}))]
    sites = tabulated(tmp_path, chr1, reads)
    assert [site.depth for site in sites[:12]] == [1] * 3 + [0] * 4 + [1] * 5
    assert not any(site.alts for site in sites)


def random_reference(rng: random.Random, length: int) -> str:
    """Bases in runs of short repeats, where indels left-align far."""
    bases = ""
    while len(bases) < length:
        unit = "".join(rng.choice("ACGT") for _ in range(rng.choice([1, 1, 2, 3])))
        bases += unit * rng.randint(1, 6)
    return bases[:length]


def random_read(rng: random.Random, name: str, chr1: str) -> AlignedSegment:
    """A read with mismatches, read Ns, indels, and skips, ending in five matching bases."""
    start = position = rng.randint(0, 30)
    bases, cigar = "", ""
    for _ in range(rng.randint(2, 7)):
        operator = rng.choices("MXnIDN", weights=[6, 2, 1, 1, 1, 1])[0]
        length = rng.randint(1, 8) if operator == "M" else rng.randint(1, 3)
        if operator in "MXn":
            for base in chr1[position : position + length]:
                bases += {
                    "M": base,
                    "n": "N",
                    "X": rng.choice([other for other in "ACGT" if other != base]),
                }[operator]
            cigar += f"{length}M"
            position += length
        elif operator == "I":
            bases += "".join(rng.choice("ACGT") for _ in range(length))
            cigar += f"{length}I"
        else:
            cigar += f"{length}{operator}"
            position += length
    bases += chr1[position : position + 5]
    return record(
        name,
        start,
        f"{cigar}5M",
        bases,
        flag=rng.choice([0, 16]),
        quals=[rng.choice([2, 20, 40, 40, 40]) for _ in bases],
        header=header_of({"chr1": chr1}),
    )


@pytest.mark.parametrize("seed", range(20))
def test_each_read_is_counted_once_where_it_holds_a_base_or_deletion(
    tmp_path: Path, seed: int
) -> None:
    rng = random.Random(seed)
    chr1 = random_reference(rng, 150)
    reads = sorted(
        (random_read(rng, f"r{index}", chr1) for index in range(30)),
        key=lambda read: read.reference_start,
    )
    for min_base_quality in (0, 30):
        sites = tabulated(tmp_path, chr1, reads, min_base_quality=min_base_quality)
        with AlignmentFile(str(tmp_path / "reads.bam")) as alignments:
            columns = list(StreamingPileupBuilder(alignments).columns("chr1", 0, len(chr1)))
        for site, spans, column in zip(sites, spanning(sites), columns, strict=True):
            assert site.depth == site.ref_reads + sum(site.alt_reads) + spans
            assert site.depth <= sum(
                entry.is_del or entry.base is not None for entry in column.pileups
            )


def test_a_read_that_spells_the_reference_across_an_insertion_and_deletion_is_a_reference_read(
    tmp_path: Path,
) -> None:
    chr1 = "ACGTAGGCTAACGTTAGCCA"
    reads = [
        record(
            "r", 0, "5M1I1D5M", chr1[0:5] + chr1[5] + chr1[6:11], header=header_of({"chr1": chr1})
        )
    ]
    sites = tabulated(tmp_path, chr1, reads)
    assert [(site.depth, site.ref_reads) for site in sites[:12]] == [(1, 1)] * 11 + [(0, 0)]
    assert not any(site.alts for site in sites)


def test_no_allele_is_left_aligned_onto_a_reference_n(tmp_path: Path) -> None:
    chr1 = "GTNAAAAGTCATG"
    reads = [record("r", 0, "6M1D5M", "GTNAAA" + "GTCAT", header=header_of({"chr1": chr1}))]
    sites = tabulated(tmp_path, chr1, reads)
    assert [site.depth for site in sites[:12]] == [1, 1] + [0] * 5 + [1] * 5
    assert not any(site.alts for site in sites)


@pytest.mark.parametrize("anchor,after,counted", [(10, 40, True), (40, 10, False)])
def test_a_deletion_is_judged_by_the_base_after_it_as_in_pileups(
    tmp_path: Path, anchor: int, after: int, counted: bool
) -> None:
    chr1 = "ACGTAGGCTAACGTTAGCCA"
    quals = [40, 40, anchor, after, 40, 40]
    reads = [
        record(
            "r", 2, "3M1D3M", chr1[2:5] + chr1[6:9], quals=quals, header=header_of({"chr1": chr1})
        )
    ]
    sites = tabulated(tmp_path, chr1, reads, min_base_quality=20)
    with AlignmentFile(str(tmp_path / "reads.bam")) as alignments:
        deleted = StreamingPileupBuilder(alignments, min_base_quality=20).pileup("chr1", 5)
    assert sites[5].depth == deleted.filtered_depth == int(counted)
    assert alleles(sites[4]) == ({"AG>A": 1} if counted else {})
