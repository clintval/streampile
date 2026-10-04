import re
from collections.abc import Iterator
from importlib.metadata import version
from pathlib import Path
from typing import Any

import pysam
import pytest
from pysam import AlignmentFile
from pysam import FastaFile

from streampile import TabulatedBase
from streampile import TabulationReader
from streampile import tabulate
from streampile._cli import main

from .records import DATA
from .records import territory
from .records import write_fasta

ROOT = Path(__file__).parent.parent

README = ROOT / "README.md"


def run(
    out: Path,
    *extra: str,
    intervals: Path = DATA / "territory.bed",
    bam: Path = DATA / "reads.bam",
    ref: Path = DATA / "reference.fa",
) -> int:
    return main([
        "tabulate",
        "--bam",
        str(bam),
        "--ref",
        str(ref),
        "--intervals",
        str(intervals),
        *extra,
        "--out",
        str(out),
    ])


def rows(path: Path) -> list[str]:
    """The rows of a plain table, without its metadata and header lines."""
    return [line for line in path.read_text().splitlines() if not line.startswith("#")]


def without_version(table: str) -> str:
    """A table without the line naming the streampile version that wrote it."""
    return re.sub(r"^##streampile-version=.*\n", "", table, flags=re.MULTILINE)


def test_tabulate_writes_the_table_shown_in_the_readme(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(ROOT)
    shown = re.search(
        r"```console\n(streampile tabulate.*?)```\n\n```text\n(.*?)```",
        README.read_text(),
        re.DOTALL,
    )
    assert shown is not None
    command = shown.group(1).replace("\\\n", " ").split()[1:]
    out = tmp_path / "counts.tsv"
    assert main([*command[: command.index("--out")], "--out", str(out)]) == 0
    written = out.read_text()
    assert f"##streampile-version={version('streampile')}\n" in written
    assert without_version(written) == without_version(shown.group(2))
    assert without_version(written) == without_version((DATA / "counts.tsv").read_text())


@pytest.mark.parametrize("name", ["counts.tsv", "counts.tsv.gz"])
def test_a_written_table_reads_back_record_for_record(tmp_path: Path, name: str) -> None:
    out = tmp_path / name
    assert run(out, "--min-base-quality", "30", "--threads", "2") == 0
    with (
        AlignmentFile(str(DATA / "reads.bam")) as reads,
        FastaFile(str(DATA / "reference.fa")) as reference,
    ):
        expected = list(
            tabulate(reads, reference, territory(("chr1", 20, 24)), min_base_quality=30)
        )
    assert list(TabulationReader.from_path(out)) == expected


def test_tabulate_reads_a_bed_territory(tmp_path: Path) -> None:
    bed = tmp_path / "territory.bed"
    bed.write_text("track name=x\n# comment\nchr1\t22\t24\tgene\nchr1\t5\t5\nchr1\t20\t22\n")
    assert run(tmp_path / "counts.tsv", intervals=bed) == 0
    assert [row.split("\t")[1] for row in rows(tmp_path / "counts.tsv")] == ["21", "22", "23", "24"]
    bed.write_text("chr1\t20\t24\nchr1\t30\t25\n")
    with pytest.raises(ValueError, match="on line 2"):
        run(tmp_path / "bad.tsv", intervals=bed)
    assert not (tmp_path / "bad.tsv").exists()


@pytest.mark.parametrize("suffix", [".gz", ".bgz", ".bgzf"])
def test_tabulate_indexes_a_bgzf_table(tmp_path: Path, suffix: str) -> None:
    out = tmp_path / f"counts.tsv{suffix}"
    assert run(out, "--index", "csi") == 0
    assert Path(f"{out}.csi").is_file()
    assert out.read_bytes()[12:16] == b"BC\x02\x00"


@pytest.mark.parametrize("index", ["tbi", "csi"])
def test_a_region_of_an_indexed_table_reads_back_as_typed_rows(tmp_path: Path, index: str) -> None:
    bed = tmp_path / "territory.bed"
    bed.write_text("chr1\t0\t60\nchr2\t0\t40\n")
    out = tmp_path / "counts.tsv.gz"
    assert run(out, "--min-base-quality", "30", "--index", index, intervals=bed) == 0
    rows = list(TabulationReader.from_path(out))
    assert list(TabulationReader.query(out, "chr1", 21, 24)) == [
        row for row in rows if row.contig == "chr1" and 22 <= row.pos <= 24
    ]
    assert list(TabulationReader.query(out, "chr2", 0, 40)) == rows[60:]
    assert any(row.alts for row in TabulationReader.query(out, "chr1", 0, 60))
    assert list(TabulationReader.query(out, "chr3", 0, 10)) == []


def test_tabulate_refuses_an_index_on_a_plain_table(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="An index and threads need a BGZF path ending in .gz,"):
        run(tmp_path / "counts.tsv", "--index", "tbi")
    assert not (tmp_path / "counts.tsv").exists()


def test_tabulate_refuses_a_sam_it_cannot_read_by_region(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="fetching by region is not available for SAM files"):
        run(tmp_path / "counts.tsv", bam=DATA / "reads.sam")
    assert list(tmp_path.iterdir()) == []


def test_tabulate_reads_a_cram_with_the_reference_it_is_given(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    with FastaFile(str(DATA / "reference.fa")) as reference:
        contigs = {name: reference.fetch(name) for name in reference.references}
    fasta = write_fasta(elsewhere / "reference.fa", contigs)
    cram = tmp_path / "reads.cram"
    with (
        AlignmentFile(str(DATA / "reads.bam")) as reads,
        AlignmentFile(str(cram), "wc", template=reads, reference_filename=str(fasta)) as sink,
    ):
        for read in reads:
            sink.write(read)
    pysam.index(str(cram))
    moved = write_fasta(tmp_path / "moved.fa", contigs)
    for path in (fasta, Path(f"{fasta}.fai")):
        path.unlink()
    monkeypatch.setenv("REF_PATH", f"{tmp_path / 'no-cache'}/%s")
    monkeypatch.setenv("REF_CACHE", f"{tmp_path / 'no-cache'}/%s")
    out = tmp_path / "counts.tsv"
    assert (
        run(out, "--min-base-quality", "30", "--min-mapping-quality", "20", bam=cram, ref=moved)
        == 0
    )
    assert rows(out) == rows(DATA / "counts.tsv")


def test_tabulate_writes_nothing_for_a_territory_it_refuses(tmp_path: Path) -> None:
    bed = tmp_path / "territory.bed"
    bed.write_text("chr1\t0\t10\nchr2\t0\t10\n")
    with FastaFile(str(DATA / "reference.fa")) as reference:
        fasta = write_fasta(tmp_path / "chr1.fa", {"chr1": reference.fetch("chr1")})
    with pytest.raises(ValueError, match="Contig chr2 is not in the reference."):
        run(tmp_path / "counts.tsv.gz", "--index", "tbi", intervals=bed, ref=fasta)
    assert not (tmp_path / "counts.tsv.gz").exists()
    assert not (tmp_path / "counts.tsv.gz.tbi").exists()


@pytest.mark.parametrize("name", ["counts.tsv", "counts.tsv.gz"])
def test_a_failed_run_leaves_no_table_and_keeps_the_old_one(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, name: str
) -> None:
    def tabulate_then_fail(*args: Any, **kwargs: Any) -> Iterator[TabulatedBase]:
        yield next(tabulate(*args, **kwargs))
        raise RuntimeError("disk full")

    monkeypatch.setattr("streampile._cli.tabulate", tabulate_then_fail)
    out = tmp_path / name
    out.write_text("old")
    index = ["--index", "tbi"] if name.endswith(".gz") else []
    with pytest.raises(RuntimeError, match="disk full"):
        run(out, *index)
    assert [path.name for path in tmp_path.iterdir()] == [name]
    assert out.read_text() == "old"


def test_tabulate_writes_to_a_path_that_is_not_a_regular_file() -> None:
    assert run(Path("/dev/null")) == 0


def test_tabulate_reads_flags_in_any_base(tmp_path: Path) -> None:
    out = tmp_path / "counts.tsv"
    assert run(out, "--exclude-flags", "0x0") == 0
    assert rows(out)[0] == "chr1\t21\tT\t5\t0\t5\t3\t2\t\t\t\t\t"


def test_a_command_is_required(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit):
        main([])
    assert "required" in capsys.readouterr().err
