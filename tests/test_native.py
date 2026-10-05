import tomllib
from collections import Counter
from importlib.metadata import version
from pathlib import Path
from typing import TYPE_CHECKING
from typing import Any

import pytest
from pysam import AlignmentFile

from streampile import Pileup
from streampile import StreamingPileupBuilder
from streampile._native import ColumnCounts
from streampile._native import sweep

from .records import DATA
from .records import write_bam
from .test_pysam_agreement import LENGTH
from .test_pysam_agreement import RANDOM_HEADER
from .test_pysam_agreement import random_reads

if TYPE_CHECKING:
    from streampile._native import NativeColumn

ROOT = Path(__file__).parent.parent

FIXTURE_SPANS = [("chr1", 0, 60), ("chr2", 0, 40)]

OPTIONS: list[dict[str, Any]] = [
    {},
    {"min_base_quality": 0},
    {"min_base_quality": 30},
    {"exclude_flags": 0},
    {"exclude_flags": 0, "min_mapping_quality": 10},
    {"proper_pairs_only": True},
    {"without_overlaps": True},
    {"without_overlaps": True, "min_base_quality": 30, "exclude_flags": 0},
]


def column(pileup: Pileup) -> "NativeColumn":
    """A Python pileup in the shape the native sweep returns."""
    entries = [
        (
            entry.alignment.query_name or "*",
            entry.alignment.flag,
            entry.pileup_type.value,
            entry.query_position,
            entry.query_position_or_next,
            entry.insertion_offset,
            entry.insertion_length,
            entry.base,
            entry.qual,
            entry.inserted_bases,
            None if (inserted := entry.inserted_qualities) is None else bytes(inserted),
        )
        for entry in pileup.pileups
    ]
    return (
        entries,
        pileup.unfiltered_depth,
        pileup.filtered_depth,
        "".join(pileup.bases),
        bytes(pileup.qualities),
    )


def python_sweep(
    path: Path, spans: list[tuple[str, int, int]], **options: Any
) -> list["NativeColumn"]:
    """Every column of the spans from the Python builder, with the native sweep's options."""
    without_overlaps = options.pop("without_overlaps", False)
    swept: list[NativeColumn] = []
    with AlignmentFile(str(path)) as reads, StreamingPileupBuilder(reads, **options) as builder:
        for contig, start, end in spans:
            for pileup in builder.columns(contig, start, end):
                swept.append(column(pileup.without_overlaps() if without_overlaps else pileup))
    return swept


@pytest.mark.parametrize("options", OPTIONS)
def test_the_rust_pileup_equals_the_python_pileup_on_the_fixture(options: dict[str, Any]) -> None:
    path = DATA / "reads.bam"
    assert sweep(path, FIXTURE_SPANS, **options) == python_sweep(path, FIXTURE_SPANS, **options)


@pytest.mark.parametrize("seed", range(25))
def test_the_rust_pileup_equals_the_python_pileup_on_random_reads(
    tmp_path: Path, seed: int
) -> None:
    reads = random_reads(seed, count=40)
    for index, read in enumerate(reads):
        read.query_name = f"t{index % 6}"
    path = write_bam(tmp_path / "reads.bam", reads, header=RANDOM_HEADER)
    spans = [("chr1", 0, LENGTH)]
    for options in OPTIONS:
        assert sweep(path, spans, **options) == python_sweep(path, spans, **options), options


def test_the_rust_sweep_refuses_what_the_python_builder_refuses(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="contig chr3 is not in the header"):
        sweep(DATA / "reads.bam", [("chr3", 0, 10)])
    with pytest.raises(ValueError, match="end 1 is before start 2"):
        sweep(DATA / "reads.bam", [("chr1", 2, 1)])
    with pytest.raises(ValueError, match="attempted to advance to chr1:0 from chr2:5"):
        sweep(DATA / "reads.bam", [("chr2", 5, 6), ("chr1", 0, 1)])
    with pytest.raises(OSError, match="missing.bam"):
        sweep(tmp_path / "missing.bam", [("chr1", 0, 1)])
    with pytest.raises(ValueError, match="contig chr3 is not in the header"):
        ColumnCounts(
            DATA / "reads.bam",
            [("chr3", 0, 1)],
            min_mapping_quality=0,
            exclude_flags=0,
            quality_floor=0,
        )


def test_the_rust_column_counts_equal_counts_of_python_pileups() -> None:
    expected: list[dict[str, int]] = []
    with AlignmentFile(str(DATA / "reads.bam")) as reads, StreamingPileupBuilder(reads) as builder:
        for contig, start, end in FIXTURE_SPANS:
            for pileup in builder.columns(contig, start, end):
                counts: Counter[str] = Counter()
                for entry in pileup.pileups:
                    if entry.is_ins:
                        counts["+"] += 1
                    elif (quality := entry.qual) is not None and quality >= 30:
                        counts["-" if entry.is_del else entry.base or "N"] += 1
                expected.append(dict(counts))
    counted = ColumnCounts(
        DATA / "reads.bam",
        FIXTURE_SPANS,
        min_mapping_quality=0,
        exclude_flags=0xF00,
        quality_floor=30,
    )
    assert list(counted) == expected


def test_the_version_comes_from_cargo() -> None:
    cargo = tomllib.loads((ROOT / "Cargo.toml").read_text())
    project = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]
    assert "version" not in project
    assert project["dynamic"] == ["version"]
    assert version("streampile") == cargo["package"]["version"]
