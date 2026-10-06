from array import array
from collections.abc import Iterator
from pathlib import Path
from types import SimpleNamespace
from typing import Any
from typing import Literal

import pytest
from pysam import AlignmentFile

from streampile import Pileup
from streampile import PileupRead
from streampile import StreamingPileupBuilder
from streampile import _native
from streampile._pileup import BASE

from .records import record
from .records import write_bam
from .test_golden import OPTIONS
from .test_golden import columns
from .test_golden import described
from .test_golden import templates
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
    """Every column of the random reads, with its templates and each entry's distances."""
    swept: list[str] = []
    with AlignmentFile(str(path)) as reads, StreamingPileupBuilder(reads, **options) as builder:
        for pileup in builder.columns("chr1", 0, LENGTH):
            swept.append(described(pileup) + templates(pileup) + ends(pileup))
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


@pytest.mark.parametrize(
    ("value", "value_type", "shown"),
    [
        (1.5, "f", "1.5"),
        (array("i", [1, 2]), None, "1,2"),
        ("M", "A", "M"),
        (4_000_000_000, "I", "4000000000"),
    ],
)
def test_a_mate_cigar_that_is_not_a_string_is_piled_up_alike_by_both_bridges(
    value: "float | str | array[int]", value_type: Literal["f", "A", "I"] | None, shown: str
) -> None:
    read = record("p", 2, "5M", "ACGTA", flag=0x1 | 0x20 | 0x40)
    read.next_reference_id = 0
    read.next_reference_start = 8
    read.set_tag("MC", value, value_type=value_type)  # pyright: ignore[reportUnknownMemberType]
    fragment = record("f", 2, "5M", "ACGTA")
    fragment.set_tag("MC", value, value_type=value_type)  # pyright: ignore[reportUnknownMemberType]
    for direct in (True, False):
        assert _native.direct_bridge(direct) is direct
        try:
            with StreamingPileupBuilder([read, fragment], exclude_flags=0) as builder:
                pileup = builder.pileup("chr1", 4)
        finally:
            assert _native.direct_bridge(True)
        assert pileup.bases == ["G", "G"]
        with pytest.raises(ValueError, match=rf"^Read p has an invalid MC tag: {shown}\.$"):
            _ = pileup.pileups[0].template_end_distance
        assert pileup.pileups[1].template_end_distance is None


def test_a_base_written_as_equals_piles_up_as_equals_on_both_bridges() -> None:
    read = record("eq", 10, "4M", "A=GT")
    for direct in (True, False):
        assert _native.direct_bridge(direct) is direct
        try:
            pileup = StreamingPileupBuilder([read]).pileup("chr1", 11)
        finally:
            assert _native.direct_bridge(True)
        assert (pileup.bases, pileup.pileups[0].base) == (["="], "=")


def test_a_record_that_does_not_fit_bam_is_refused_naming_the_read() -> None:
    read = record("big", 10, "4M", "ACGT")
    attributes = {
        name: getattr(read, name)
        for name in (
            "query_name",
            "cigartuples",
            "query_sequence",
            "query_qualities",
            "reference_id",
            "mapping_quality",
            "flag",
            "next_reference_id",
            "next_reference_start",
            "template_length",
        )
    }

    def has_tag(_tag: str) -> bool:
        return False

    distant = SimpleNamespace(**attributes, reference_start=2**40, has_tag=has_tag)
    with pytest.raises(ValueError, match=r"^Read big is invalid: its position is too long\.$"):
        PileupRead._make([distant, 0, 0, BASE])
