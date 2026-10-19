import random
from collections.abc import Iterable
from io import StringIO
from pathlib import Path
from typing import cast
from typing import final

import pytest
from pysam import AlignedSegment
from pysam import AlignmentFile

from streampile import DEFAULT_EXCLUDE_FLAGS
from streampile import ContigBins
from streampile import _native
from streampile import bin_templates
from streampile import write_wig

from .records import DATA
from .records import HEADER
from .records import record
from .records import write_bam


def binned(
    tmp_path: Path,
    reads: Iterable[AlignedSegment],
    window: int = 100,
    *,
    min_mapping_quality: int = 0,
    exclude_flags: int = DEFAULT_EXCLUDE_FLAGS,
) -> dict[str, tuple[int, ...]]:
    """The templates in each window of each contig, from reads written to an indexed BAM."""
    path = write_bam(tmp_path / "reads.bam", reads)
    with AlignmentFile(str(path)) as alignments:
        bins = bin_templates(
            alignments,
            window,
            min_mapping_quality=min_mapping_quality,
            exclude_flags=exclude_flags,
        )
        return {contig.contig: contig.templates for contig in bins}


def mate(
    name: str,
    start: int,
    flag: int,
    mate_start: int,
    *,
    contig: str = "chr1",
    mate_contig: str = "chr1",
    mapq: int = 60,
) -> AlignedSegment:
    """One read of a pair, whose mate starts at `mate_start` of `mate_contig`."""
    read = record(name, start, "10M", "ACGTACGTAC", contig=contig, flag=flag, mapq=mapq)
    read.next_reference_name = mate_contig
    read.next_reference_start = mate_start
    return read


def test_bin_templates_counts_the_fixture_by_hand() -> None:
    with AlignmentFile(str(DATA / "reads.bam")) as reads:
        assert list(bin_templates(reads, 20)) == [
            ContigBins(contig="chr1", window=20, templates=(5, 3, 1)),
            ContigBins(contig="chr2", window=20, templates=(1, 0)),
        ]
    with AlignmentFile(str(DATA / "reads.bam")) as reads:
        counted = [contig.templates for contig in bin_templates(reads, 20, exclude_flags=0x900)]
    assert counted == [(5, 4, 1), (1, 0)]


def test_the_two_reads_of_a_pair_count_once_in_the_window_of_the_first(tmp_path: Path) -> None:
    reads = [
        mate("apart", 95, 99, 105),
        mate("apart", 105, 147, 95),
        mate("level", 150, 99, 150),
        mate("level", 150, 147, 150),
    ]
    assert binned(tmp_path, reads)["chr1"][:2] == (1, 1)


def test_a_pair_with_one_read_filtered_out_counts_once_at_the_other(tmp_path: Path) -> None:
    reads = [
        mate("first", 50, 99, 150, mapq=5),
        mate("second", 60, 99, 260),
        mate("first", 150, 147, 50),
        mate("second", 260, 147, 60, mapq=5),
    ]
    assert binned(tmp_path, reads, min_mapping_quality=20)["chr1"][:3] == (1, 1, 0)


def test_a_pair_across_contigs_counts_once_on_the_contig_binned_first(tmp_path: Path) -> None:
    reads = [
        mate("split", 100, 65, 50, mate_contig="chr2"),
        mate("split", 50, 129, 100, contig="chr2"),
    ]
    path = write_bam(tmp_path / "reads.bam", reads)

    def counted(contigs: list[str] | None) -> dict[str, int]:
        with AlignmentFile(str(path)) as alignments:
            bins = bin_templates(alignments, 1000, contigs=contigs)
            return {contig.contig: contig.templates[0] for contig in bins}

    assert counted(None) == {"chr1": 1, "chr2": 0}
    assert counted(["chr2", "chr1"]) == {"chr2": 1, "chr1": 0}
    assert counted(["chr2"]) == {"chr2": 1}
    assert counted(["chr1"]) == {"chr1": 1}


def test_a_read_whose_mate_is_unmapped_counts_alone(tmp_path: Path) -> None:
    placed = record("orphan", 100, "*", "ACGT", flag=133)
    placed.next_reference_name = "chr1"
    placed.next_reference_start = 100
    reads = [mate("orphan", 100, 73, 100), placed]
    assert binned(tmp_path, reads)["chr1"][1] == 1


def test_secondary_and_supplementary_records_never_count(tmp_path: Path) -> None:
    reads = [
        record("chimera", 10, "10M", "ACGTACGTAC"),
        record("chimera", 30, "10M", "ACGTACGTAC", flag=256),
        record("chimera", 50, "10M", "ACGTACGTAC", flag=2048),
        record("copy", 70, "10M", "ACGTACGTAC", flag=1024),
    ]
    assert binned(tmp_path, reads, exclude_flags=0x900)["chr1"][0] == 2


@pytest.mark.parametrize("flags", [0, 0x100, 0x800, 0x400])
def test_exclude_flags_must_hold_the_secondary_and_supplementary_flags(flags: int) -> None:
    message = (
        f"exclude_flags must hold the secondary and supplementary flags, 0x900, found: {flags:#x}"
    )
    with (
        AlignmentFile(str(DATA / "reads.bam")) as reads,
        pytest.raises(ValueError, match=rf"^{message}$"),
    ):
        bin_templates(reads, 20, exclude_flags=flags)


def test_windows_are_half_open_and_the_last_is_shorter(tmp_path: Path) -> None:
    reads = [record(name, start, "1M", "A") for name, start in (("a", 299), ("b", 300), ("c", 999))]
    assert binned(tmp_path, reads, window=300) == {"chr1": (1, 1, 0, 1), "chr2": (0, 0, 0, 0)}


def test_a_window_longer_than_a_contig_is_the_contig() -> None:
    with AlignmentFile(str(DATA / "reads.bam")) as reads:
        bins = list(bin_templates(reads, 50))
    assert [(contig.window, contig.templates) for contig in bins] == [(50, (9, 0)), (40, (1,))]
    wig = StringIO()
    write_wig(bins, wig)
    assert wig.getvalue() == (
        "fixedStep chrom=chr1 start=1 step=50 span=50\n9\n0\n"
        + "fixedStep chrom=chr2 start=1 step=40 span=40\n1\n"
    )


def test_bin_templates_floors_mapping_quality(tmp_path: Path) -> None:
    reads = [record(name, 10, "4M", "ACGT", mapq=mapq) for name, mapq in (("a", 19), ("b", 20))]
    reads.append(record("c", 10, "4M", "ACGT", mapq=255))
    assert binned(tmp_path, reads, min_mapping_quality=20)["chr1"][0] == 2


def test_bin_templates_checks_contigs_when_called() -> None:
    with AlignmentFile(str(DATA / "reads.bam")) as reads:
        with pytest.raises(ValueError, match="^Contig chr3 is not in the alignment header.$"):
            bin_templates(reads, 20, contigs=["chr1", "chr3"])
        with pytest.raises(ValueError, match="^Contig chr1 is given more than once.$"):
            bin_templates(reads, 20, contigs=["chr1", "chr2", "chr1"])


@pytest.mark.parametrize(
    ("window", "min_mapping_quality", "exclude_flags", "message"),
    [
        (0, 0, 0xF00, "window must be from 1 to 4294967295, found: 0"),
        (2**32, 0, 0xF00, "window must be from 1 to 4294967295, found: 4294967296"),
        (20, 256, 0xF00, "min_mapping_quality must be from 0 to 255, found: 256"),
        (20, 0, 0x10000, "exclude_flags must be from 0 to 65535, found: 65536"),
    ],
)
def test_bin_templates_refuses_options_out_of_range(
    window: int, min_mapping_quality: int, exclude_flags: int, message: str
) -> None:
    with (
        AlignmentFile(str(DATA / "reads.bam")) as reads,
        pytest.raises(ValueError, match=rf"^{message}$"),
    ):
        bin_templates(
            reads, window, min_mapping_quality=min_mapping_quality, exclude_flags=exclude_flags
        )


@final
class Unsorted:
    """Alignments that give a contig's reads as listed, in or out of order."""

    def __init__(self, reads: list[AlignedSegment]) -> None:
        self.reads = reads

    def fetch(self, contig: str) -> list[AlignedSegment]:
        return [read for read in self.reads if read.reference_name == contig]


def unsorted(reads: list[AlignedSegment]) -> AlignmentFile:
    """Reads listed in any order, as the alignments a `TemplateBins` fetches from."""
    return cast(AlignmentFile, cast(object, Unsorted(reads)))


def test_bin_templates_refuses_reads_out_of_order_or_past_the_end() -> None:
    reads = [record("late", 50, "4M", "ACGT"), record("early", 10, "4M", "ACGT")]
    with pytest.raises(ValueError, match="^Records are out of coordinate order at early.$"):
        list(_native.TemplateBins(unsorted(reads), [("chr1", 0, 1000)], 100))
    reads = [record("over", 500, "4M", "ACGT")]
    with pytest.raises(
        ValueError, match="^Read over starts past the end of chr1, which has 400 bases.$"
    ):
        list(_native.TemplateBins(unsorted(reads), [("chr1", 0, 400)], 100))


def oracle(
    path: Path, window: int, contigs: list[str], min_mapping_quality: int
) -> list[list[int]]:
    """The templates in each window, counting each read name once, at its first passing read."""
    counted: set[str] = set()
    bins: list[list[int]] = []
    with AlignmentFile(str(path)) as alignments:
        for contig in contigs:
            windows = [0] * -(-alignments.get_reference_length(contig) // window)
            for read in alignments.fetch(contig):
                passes = not read.is_unmapped and (read.flag & DEFAULT_EXCLUDE_FLAGS) == 0
                if not passes or read.mapping_quality < min_mapping_quality:
                    continue
                if read.query_name not in counted:
                    counted.add(read.query_name or "")
                    windows[read.reference_start // window] += 1
            bins.append(windows)
    return bins


@pytest.mark.parametrize("seed", range(20))
def test_bin_templates_counts_each_name_once_on_random_reads(tmp_path: Path, seed: int) -> None:
    rng = random.Random(seed)
    contigs = ["chr1", "chr2"]
    reads: list[AlignedSegment] = []
    for index in range(300):
        name = f"t{index}"
        contig = rng.choice(contigs)
        start = rng.randrange(990)
        if rng.random() < 0.3:
            reads.append(
                record(name, start, "10M", "ACGTACGTAC", contig=contig, mapq=rng.randrange(61))
            )
            continue
        other = contig if rng.random() < 0.8 else rng.choice(contigs)
        other_start = rng.randrange(990)
        first, second = (65, 129) if rng.random() < 0.5 else (129, 65)
        flags = (first, second) if rng.random() < 0.9 else (first | 0x400, second | 0x400)
        reads.append(
            mate(
                name,
                start,
                flags[0],
                other_start,
                contig=contig,
                mate_contig=other,
                mapq=rng.randrange(61),
            )
        )
        reads.append(
            mate(
                name,
                other_start,
                flags[1],
                start,
                contig=other,
                mate_contig=contig,
                mapq=rng.randrange(61),
            )
        )
    reads.sort(key=lambda read: (read.reference_id, read.reference_start))
    path = write_bam(tmp_path / "reads.bam", reads, HEADER)
    order = contigs if rng.random() < 0.5 else contigs[::-1]
    window = rng.choice([1, 7, 100, 1000, 5000])
    with AlignmentFile(str(path)) as alignments:
        ours = [
            list(contig.templates)
            for contig in bin_templates(alignments, window, contigs=order, min_mapping_quality=20)
        ]
    assert ours == oracle(path, window, order, 20)
