import gzip
from pathlib import Path

import pybgzf
import pytest
from pybgzf import IndexFormat

from streampile import TabulatedBase
from streampile import TabulationReader
from streampile import TabulationWriter

BASES = [
    TabulatedBase(contig="chr1", pos=1, ref="A", depth=0, ref_reads=0),
    TabulatedBase(
        contig="chr1",
        pos=2,
        ref="C",
        depth=9,
        ref_reads=5,
        alt_refs=("C", "CA"),
        alts=("T", "C"),
        alt_reads=(3, 1),
    ),
    TabulatedBase(contig="chr1", pos=3, ref="A", depth=9, ref_reads=8),
    TabulatedBase(
        contig="chr2",
        pos=7,
        ref="G",
        depth=4,
        ref_reads=3,
        alt_refs=("G",),
        alts=("GTT",),
        alt_reads=(1,),
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
    assert lines[0] == "contig\tpos\tref\tdepth\tref_reads\talt_refs\talts\talt_reads"
    assert lines[1] == "chr1\t1\tA\t0\t0\t\t\t"
    assert lines[2] == "chr1\t2\tC\t9\t5\tC,CA\tT,C\t3,1"
    assert list(TabulationReader.from_path(tmp_path / "bases.tsv")) == BASES


@pytest.mark.parametrize("suffix", [".gz", ".bgz"])
def test_a_compressed_table_is_bgzf_with_an_end_of_file_block(tmp_path: Path, suffix: str) -> None:
    path = tmp_path / f"bases.tsv{suffix}"
    write(path, BASES, threads=2)
    raw = path.read_bytes()
    assert raw[12:16] == b"BC\x02\x00"
    assert raw.endswith(bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000"))
    with gzip.open(path, "rt") as handle:
        assert handle.readline().startswith("contig\tpos")
    assert list(TabulationReader.from_path(path)) == BASES


@pytest.mark.parametrize("index", [IndexFormat.TBI, IndexFormat.CSI])
def test_an_index_finds_the_rows_of_a_region(tmp_path: Path, index: IndexFormat) -> None:
    path = tmp_path / "bases.tsv.gz"
    write(path, BASES, index=index)
    with pybgzf.IndexedReader(path) as reader:
        assert list(reader.query("chr1", 1, 3)) == [
            "chr1\t2\tC\t9\t5\tC,CA\tT,C\t3,1",
            "chr1\t3\tA\t9\t8\t\t\t",
        ]
        assert list(reader.query("chr2", 0, 100)) == ["chr2\t7\tG\t4\t3\tG\tGTT\t1"]


def test_an_indexed_table_refuses_rows_out_of_order(tmp_path: Path) -> None:
    with pytest.raises(Exception, match="(?i)sort|order"):
        write(tmp_path / "bases.tsv.gz", [BASES[2], BASES[0]], index=IndexFormat.TBI)


def test_the_writer_refuses_bad_options_before_opening_a_file(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="An index and threads need a BGZF path"):
        TabulationWriter.from_path(tmp_path / "bases.tsv", threads=2)
    with pytest.raises(ValueError, match="An index_path needs an index"):
        TabulationWriter.from_path(tmp_path / "bases.tsv.gz", index_path=tmp_path / "bases.csi")
    with pytest.raises(ValueError):
        TabulationWriter.from_path(tmp_path / "bases.tsv.gz", rename={"nothing": "x"})
    assert list(tmp_path.iterdir()) == []
