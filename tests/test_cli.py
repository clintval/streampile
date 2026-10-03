import re
from pathlib import Path

import pytest
from pysam import AlignmentFile
from pysam import FastaFile

from streampile import TabulationReader
from streampile import tabulate
from streampile._cli import main

from .records import DATA

README = Path(__file__).parent.parent / "README.md"


def run(out: Path, *extra: str) -> int:
    return main([
        "tabulate",
        "--bam",
        str(DATA / "reads.bam"),
        "--ref",
        str(DATA / "reference.fa"),
        "--intervals",
        str(DATA / "territory.bed"),
        *extra,
        "--out",
        str(out),
    ])


def test_tabulate_writes_the_table_shown_in_the_readme(tmp_path: Path) -> None:
    out = tmp_path / "counts.tsv"
    assert run(out, "--min-base-quality", "30", "--min-mapping-quality", "20") == 0
    shown = re.search(r"```text\n(contig.*?)```", README.read_text(), re.DOTALL)
    assert shown is not None
    assert out.read_text() == shown.group(1) == (DATA / "counts.tsv").read_text()


@pytest.mark.parametrize("name", ["counts.tsv", "counts.tsv.gz"])
def test_a_written_table_reads_back_record_for_record(tmp_path: Path, name: str) -> None:
    out = tmp_path / name
    assert run(out, "--min-base-quality", "30", "--threads", "2") == 0
    with (
        AlignmentFile(str(DATA / "reads.bam")) as reads,
        FastaFile(str(DATA / "reference.fa")) as reference,
    ):
        expected = list(tabulate(reads, reference, [("chr1", 20, 24)], min_base_quality=30))
    assert list(TabulationReader.from_path(out)) == expected


def test_tabulate_indexes_a_bgzf_table(tmp_path: Path) -> None:
    assert run(tmp_path / "counts.tsv.gz", "--index", "csi") == 0
    assert (tmp_path / "counts.tsv.gz.csi").is_file()


def test_tabulate_refuses_an_index_on_a_plain_table(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="An index and threads need a BGZF path"):
        run(tmp_path / "counts.tsv", "--index", "tbi")
    assert not (tmp_path / "counts.tsv").exists()


def test_tabulate_reads_flags_in_any_base(tmp_path: Path) -> None:
    out = tmp_path / "counts.tsv"
    assert run(out, "--exclude-flags", "0x0") == 0
    assert out.read_text().splitlines()[1] == "chr1\t21\tT\t5\t5\t\t\t"


def test_a_command_is_required(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit):
        main([])
    assert "required" in capsys.readouterr().err
