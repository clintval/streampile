from pathlib import Path
from typing import Any

import pytest
from pysam import AlignmentFile

from streampile import Pileup
from streampile import StreamingPileupBuilder
from streampile._cli import main

from .records import DATA

GOLDEN = DATA / "golden"

OPTIONS: dict[str, dict[str, Any]] = {
    "default": {},
    "everything": {"exclude_flags": 0, "min_base_quality": 0},
    "q30": {"min_base_quality": 30},
    "proper_pairs_mq10": {"proper_pairs_only": True, "min_mapping_quality": 10},
}

TABLES: dict[str, list[str]] = {
    "default": [],
    "q30_mq20": ["--min-base-quality", "30", "--min-mapping-quality", "20"],
    "q13_everything": ["--min-base-quality", "13", "--exclude-flags", "0"],
}


def described(pileup: Pileup) -> str:
    """A pileup's depths, views, and every field of every entry, on one line."""
    entries = [
        "/".join(
            str(field)
            for field in (
                entry.alignment.query_name,
                entry.alignment.flag,
                entry.pileup_type.value,
                entry.query_position,
                entry.query_position_or_next,
                entry.insertion_offset,
                entry.insertion_length,
                entry.base,
                entry.qual,
                entry.inserted_bases,
                entry.inserted_qualities,
            )
        ).replace(" ", "")
        for entry in pileup.pileups
    ]
    return "\t".join([
        str(pileup.unfiltered_depth),
        str(pileup.filtered_depth),
        "".join(pileup.bases),
        ",".join(map(str, pileup.qualities)),
        " ".join(entries),
    ])


def columns(options: dict[str, Any]) -> str:
    lines: list[str] = []
    with (
        AlignmentFile(str(DATA / "reads.bam")) as reads,
        StreamingPileupBuilder(reads, **options) as builder,
    ):
        for contig, end in (("chr1", 60), ("chr2", 40)):
            for pileup in builder.columns(contig, 0, end):
                position = f"{pileup.reference_name}\t{pileup.reference_pos}"
                lines.append(f"{position}\tall\t{described(pileup)}")
                lines.append(f"{position}\tkept\t{described(pileup.without_overlaps())}")
    return "\n".join(lines) + "\n"


def rows(table: str) -> list[str]:
    """The header and rows of a table, without the metadata that names its inputs and version."""
    return [line for line in table.splitlines() if not line.startswith("##")]


@pytest.mark.parametrize("name", OPTIONS)
def test_the_columns_of_the_fixture_match_their_golden_file(name: str) -> None:
    assert columns(OPTIONS[name]) == (GOLDEN / f"columns_{name}.tsv").read_text()


@pytest.mark.parametrize("name", TABLES)
def test_the_table_of_the_fixture_matches_its_golden_file(tmp_path: Path, name: str) -> None:
    bed = tmp_path / "territory.bed"
    bed.write_text("chr1\t0\t60\nchr2\t0\t40\n")
    out = tmp_path / "table.tsv"
    command = ["tabulate", "--bam", str(DATA / "reads.bam"), "--ref", str(DATA / "reference.fa")]
    assert main([*command, "--intervals", str(bed), *TABLES[name], "--out", str(out)]) == 0
    assert rows(out.read_text()) == rows((GOLDEN / f"table_{name}.tsv").read_text())
