import gzip
from dataclasses import replace
from importlib.metadata import version
from pathlib import Path

import pybgzf
import pytest
from pybgzf import IndexFormat
from typeline import Comment

from streampile import TabulatedBase
from streampile import TabulationReader
from streampile import TabulationWriter

VERSION = version("streampile")

HEADER = "\t".join([
    "#contig",
    "pos",
    "ref",
    "depth",
    "ref_reads",
    "ref_fwd",
    "ref_rev",
    "alt_refs",
    "alts",
    "alt_reads",
    "alt_fwd",
    "alt_rev",
])

BASES = [
    TabulatedBase(contig="chr1", pos=1, ref="A", depth=0, ref_reads=0, ref_fwd=0, ref_rev=0),
    TabulatedBase(
        contig="chr1",
        pos=2,
        ref="C",
        depth=9,
        ref_reads=5,
        ref_fwd=3,
        ref_rev=2,
        alt_refs=("C", "CA"),
        alts=("T", "C"),
        alt_reads=(3, 1),
        alt_fwd=(2, 0),
        alt_rev=(1, 1),
    ),
    TabulatedBase(contig="chr1", pos=3, ref="A", depth=9, ref_reads=8, ref_fwd=4, ref_rev=4),
    TabulatedBase(
        contig="chr2",
        pos=7,
        ref="G",
        depth=4,
        ref_reads=3,
        ref_fwd=3,
        ref_rev=0,
        alt_refs=("G",),
        alts=("GTT",),
        alt_reads=(1,),
        alt_fwd=(0,),
        alt_rev=(1,),
    ),
]


def write(
    path: Path, bases: list[TabulatedBase], index: IndexFormat | None = None, threads: int = 1
) -> None:
    with TabulationWriter.from_path(path, index=index, threads=threads) as writer:
        writer.write_header()
        for base in bases:
            writer.write(base)


def test_a_table_round_trips_with_empty_alleles(tmp_path: Path) -> None:
    write(tmp_path / "bases.tsv", BASES)
    lines = (tmp_path / "bases.tsv").read_text().splitlines()
    assert lines[:2] == ["##streampile-tabulation=1", f"##streampile-version={VERSION}"]
    assert lines[2] == HEADER
    assert lines[3] == "chr1\t1\tA\t0\t0\t0\t0\t\t\t\t\t"
    assert lines[4] == "chr1\t2\tC\t9\t5\t3\t2\tC,CA\tT,C\t3,1\t2,0\t1,1"
    assert list(TabulationReader.from_path(tmp_path / "bases.tsv")) == BASES


@pytest.mark.parametrize("suffix", [".gz", ".bgz"])
def test_a_compressed_table_is_bgzf_with_an_end_of_file_block(tmp_path: Path, suffix: str) -> None:
    path = tmp_path / f"bases.tsv{suffix}"
    write(path, BASES, threads=2)
    raw = path.read_bytes()
    assert raw[12:16] == b"BC\x02\x00"
    assert raw.endswith(bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000"))
    with gzip.open(path, "rt") as handle:
        assert handle.readline() == "##streampile-tabulation=1\n"
    assert list(TabulationReader.from_path(path)) == BASES


@pytest.mark.parametrize("index", [IndexFormat.TBI, IndexFormat.CSI])
def test_an_index_finds_the_rows_of_a_region(tmp_path: Path, index: IndexFormat) -> None:
    path = tmp_path / "bases.tsv.gz"
    write(path, BASES, index=index)
    with pybgzf.IndexedReader(path) as reader:
        assert list(reader.query("chr1", 1, 3)) == [
            "chr1\t2\tC\t9\t5\t3\t2\tC,CA\tT,C\t3,1\t2,0\t1,1",
            "chr1\t3\tA\t9\t8\t4\t4\t\t\t\t\t",
        ]
        assert list(reader.query("chr2", 0, 100)) == ["chr2\t7\tG\t4\t3\t3\t0\tG\tGTT\t1\t0\t1"]


def test_an_indexed_table_refuses_rows_out_of_order(tmp_path: Path) -> None:
    with pytest.raises(Exception, match="(?i)sort|order"):
        write(tmp_path / "bases.tsv.gz", [BASES[2], BASES[0]], index=IndexFormat.TBI)


@pytest.mark.parametrize("index", [None, IndexFormat.TBI, IndexFormat.CSI])
def test_the_metadata_and_header_come_first_and_once(
    tmp_path: Path, index: IndexFormat | None
) -> None:
    path = tmp_path / "bases.tsv.gz"
    metadata = {"min_base_quality": 30, "territory": "territory.bed"}
    with TabulationWriter.from_path(path, index=index, metadata=metadata) as writer:
        writer.write_comment("made by streampile")
        writer.write_header()
        writer.write(BASES[0])
        writer.write_header()
        for base in BASES[1:]:
            writer.write(base)
    with gzip.open(path, "rt") as handle:
        lines = handle.read().splitlines()
    assert lines[:6] == [
        "##streampile-tabulation=1",
        f"##streampile-version={VERSION}",
        "##min_base_quality=30",
        "##territory=territory.bed",
        HEADER,
        "## made by streampile",
    ]
    assert len(lines) == 6 + len(BASES)
    if index is not None:
        with pybgzf.IndexedReader(path) as reader:
            assert list(reader.query("chr1", 0, 1)) == ["chr1\t1\tA\t0\t0\t0\t0\t\t\t\t\t"]
    comments: list[Comment] = []
    table = TabulationReader.from_path(path, on_comment=comments.append)
    assert table.metadata == {
        "streampile-tabulation": "1",
        "streampile-version": VERSION,
        "min_base_quality": "30",
        "territory": "territory.bed",
    }
    assert list(table) == BASES
    assert [comment.text for comment in comments][-1] == "## made by streampile"


def test_a_reader_keeps_the_columns_a_later_version_appends(tmp_path: Path) -> None:
    path = tmp_path / "bases.tsv"
    write(path, BASES[:2])
    lines = path.read_text().splitlines()
    lines[2] += "\tlater"
    lines[3] += "\t7"
    lines[4] += "\t8"
    path.write_text("\n".join(lines) + "\n")
    later = list(TabulationReader.from_path(path))
    assert later == [replace(BASES[0], extra=("7",)), replace(BASES[1], extra=("8",))]
    write(tmp_path / "again.tsv", later)
    assert list(TabulationReader.from_path(tmp_path / "again.tsv")) == later


def test_a_reader_refuses_another_format_version(tmp_path: Path) -> None:
    path = tmp_path / "bases.tsv"
    write(path, BASES)
    path.write_text(
        path.read_text().replace("##streampile-tabulation=1", "##streampile-tabulation=2")
    )
    with pytest.raises(ValueError, match="in format version 2, but this reader reads version 1"):
        TabulationReader.from_path(path)


def test_the_writer_refuses_bad_options_before_opening_a_file(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="An index and threads need a BGZF path"):
        TabulationWriter.from_path(tmp_path / "bases.tsv", threads=2)
    with pytest.raises(ValueError, match="An index_path needs an index"):
        TabulationWriter.from_path(tmp_path / "bases.tsv.gz", index_path=tmp_path / "bases.csi")
    with pytest.raises(ValueError):
        TabulationWriter.from_path(tmp_path / "bases.tsv.gz", rename={"nothing": "x"})
    assert list(tmp_path.iterdir()) == []
