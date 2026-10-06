from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pytest
from pysam import AlignmentFile

from streampile import Pileup
from streampile import StreamingPileupBuilder
from streampile import _native

from .records import write_bam
from .test_golden import OPTIONS
from .test_golden import columns
from .test_golden import described
from .test_pysam_agreement import LENGTH
from .test_pysam_agreement import RANDOM_HEADER
from .test_pysam_agreement import random_reads


@pytest.fixture
def through_attributes() -> Iterator[None]:
    assert not _native.direct_bridge(False)
    try:
        yield
    finally:
        assert _native.direct_bridge(True)


def test_records_are_read_from_htslib_directly_on_the_pinned_pysam() -> None:
    assert _native.direct_bridge()


@pytest.mark.usefixtures("through_attributes")
@pytest.mark.parametrize("name", OPTIONS)
def test_records_read_through_their_attributes_pile_up_the_same(name: str) -> None:
    assert not _native.direct_bridge()
    assert (
        columns(OPTIONS[name])
        == (Path(__file__).parent / f"data/golden/columns_{name}.tsv").read_text()
    )


def sweep(path: Path, **options: Any) -> list[str]:
    """Every column of the random reads, with each entry's distances to both fragment ends."""
    swept: list[str] = []
    with AlignmentFile(str(path)) as reads, StreamingPileupBuilder(reads, **options) as builder:
        for pileup in builder.columns("chr1", 0, LENGTH):
            swept.append(described(pileup) + described(pileup.without_overlaps()) + ends(pileup))
    return swept


def ends(pileup: Pileup) -> str:
    distances: list[object] = []
    for entry in pileup.pileups:
        try:
            distances.append((entry.five_prime_distance, entry.template_end_distance))
        except ValueError as error:
            distances.append(str(error))
    return repr(distances)


@pytest.mark.parametrize("seed", range(10))
def test_both_bridges_read_random_records_alike(tmp_path: Path, seed: int) -> None:
    reads = random_reads(seed, count=40)
    for index, read in enumerate(reads):
        read.query_name = f"t{index % 6}"
        if read.is_paired and index % 3:
            read.next_reference_id = 0
            read.next_reference_start = read.reference_start + 5
            read.set_tag("MC", "10M" if index % 3 == 1 else "2Q")  # pyright: ignore[reportUnknownMemberType]
    path = write_bam(tmp_path / "reads.bam", reads, header=RANDOM_HEADER)
    direct = sweep(path, exclude_flags=0, min_base_quality=20)
    assert not _native.direct_bridge(False)
    try:
        assert sweep(path, exclude_flags=0, min_base_quality=20) == direct
    finally:
        assert _native.direct_bridge(True)
