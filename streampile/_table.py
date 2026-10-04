from collections.abc import Iterator
from collections.abc import Mapping
from dataclasses import dataclass
from importlib.metadata import version
from inspect import unwrap
from io import StringIO
from pathlib import Path
from types import MappingProxyType
from typing import Final
from typing import TextIO

import pybgzf
from pybgzf import BGZF_SUFFIXES
from pybgzf import Columns
from pybgzf import IndexFormat
from typeline import Codecs
from typeline import Comment
from typeline import ExtraColumns
from typeline import FixedRecordType
from typeline import ReaderOptions
from typeline import SubscriptableClassmethod
from typeline import TsvReader
from typeline import TsvWriter
from typeline import WriterOptions
from typeline.codecs import delimited
from typing_extensions import Self
from typing_extensions import Unpack
from typing_extensions import override


@dataclass(frozen=True, kw_only=True)
class TabulatedBase:
    """The reads at one reference base, one row of a tabulation.

    A read is informative at a base, and counted in `depth`, when it holds an aligned base there
    at the base-quality floor, or when the base lies inside an allele the read was counted for.
    Each informative read is counted once: as a reference read, as a read of the one allele it
    has anchored at this base, or, when an allele anchored at an earlier base spans this one,
    in `depth` alone. So `depth` is `ref_reads`, plus `alt_reads`, plus the reads with an allele
    spanning the base. A read that skips over the base, with the CIGAR `N` operator, observed no
    base there, so it is neither a reference read nor a read of an allele, and is not in `depth`,
    unlike in a pileup's `unfiltered_depth`. A read holding an `N` base there, a no-call, is not
    informative either, and is counted in `no_calls` instead, so the molecular depth at the base
    is `depth + no_calls`, less any read left out by the base-quality floor or for an allele it
    could not be counted for.

    The alleles anchored at the base are normalized VCF alleles at `pos`, held in parallel
    tuples, most reads first: allele `i` has reference bases `alt_refs[i]`, alternate bases
    `alts[i]`, and `alt_reads[i]` reads, e.g. `C`, `T` for an SNV, `CA`, `C` for a deletion, or
    `C`, `CT` for an insertion.

    Reads are split by the strand they are mapped to, so `ref_reads` is `ref_fwd + ref_rev`, and
    `alt_reads[i]` is `alt_fwd[i] + alt_rev[i]`.

    Within a format version, columns are only ever appended, each keeping its meaning, so a table
    written by a later streampile reads here, with the columns this version does not know kept,
    as text, in `extra`.

    Attributes:
        contig: the name of the contig.
        pos: the 1-based position of the base, as in VCF.
        ref: the upper-cased reference base.
        depth: the number of informative reads.
        no_calls: the number of reads holding an `N` base at the base, which are not in `depth`.
        ref_reads: the number of informative reads with no allele at or spanning the base.
        ref_fwd: the reference reads mapped to the forward strand.
        ref_rev: the reference reads mapped to the reverse strand.
        alt_refs: the reference bases of each allele.
        alts: the alternate bases of each allele.
        alt_reads: the number of reads of each allele.
        alt_fwd: the reads of each allele mapped to the forward strand.
        alt_rev: the reads of each allele mapped to the reverse strand.
        extra: the columns after the known ones, as text, from a table written by a later
            version within the same format version, written back after the known columns.
    """

    contig: str
    pos: int
    ref: str
    depth: int
    no_calls: int
    ref_reads: int
    ref_fwd: int
    ref_rev: int
    alt_refs: tuple[str, ...] = ()
    alts: tuple[str, ...] = ()
    alt_reads: tuple[int, ...] = ()
    alt_fwd: tuple[int, ...] = ()
    alt_rev: tuple[int, ...] = ()
    extra: ExtraColumns = ()


TABULATION_CODECS: Final[Codecs] = {
    tuple[str, ...]: delimited(str, container=tuple),
    tuple[int, ...]: delimited(int, container=tuple),
}
"""How the allele fields of a `TabulatedBase` are read and written, comma-separated."""

TABULATION_FORMAT: Final[str] = "1"
"""The version of the table's format, within which columns are only ever appended."""

TABULATION_COLUMNS: Final[Columns] = Columns(
    refname=1, start=2, end=None, zero_based=False, meta_char="#"
)
"""Where an index finds each row: the contig and 1-based position, past the `#` lines."""

TABULATION_RENAME: Final[Mapping[str, str]] = MappingProxyType({"contig": "#contig"})
"""The header names its first column `#contig`, so an index and a VCF-minded reader skip it."""

METADATA_PREFIX: Final[str] = "##"
"""The prefix of the metadata lines before the header, `##key=value`."""

STREAMPILE_VERSION: Final[str] = version("streampile")
"""The version of streampile, written in the metadata of every table."""


class TabulationReader(TsvReader[TabulatedBase], FixedRecordType):
    """A reader of tabulated bases, from a TSV with a header, gzipped or BGZF or not.

    Attributes:
        metadata: the `##key=value` lines before the header, such as the table's format version
            under `streampile-tabulation` and the parameters it was tabulated with.
    """

    @override
    def __init__(self, handle: TextIO, /, **options: Unpack[ReaderOptions]) -> None:
        """Start a reader with the tabulation's codecs and header unless others are given.

        Args:
            handle: a file-like object to read the table from.
            options: the options of the reader, with tabulation defaults for any not given.

        Raises:
            ValueError: if the table declares a format version other than this reader's.
        """
        _ = options.setdefault("codecs", TABULATION_CODECS)
        _ = options.setdefault("quoting", False)
        _ = options.setdefault("rename", TABULATION_RENAME)
        _ = options.setdefault("comment_prefixes", (METADATA_PREFIX,))
        self.metadata: dict[str, str] = {}
        forward = options.get("on_comment")
        reading_metadata = True

        def on_comment(comment: Comment) -> None:
            key, found, value = comment.text.removeprefix(METADATA_PREFIX).partition("=")
            if reading_metadata and found:
                self.metadata[key] = value
            if forward is not None:
                forward(comment)

        options["on_comment"] = on_comment
        super().__init__(handle, **options)
        reading_metadata = False
        found_format = self.metadata.get("streampile-tabulation", TABULATION_FORMAT)
        if found_format != TABULATION_FORMAT:
            self.close()
            raise ValueError(
                f"The table is in format version {found_format}, but this reader reads version"
                + f" {TABULATION_FORMAT}."
            )

    @classmethod
    def query(cls, path: Path | str, refname: str, start: int, end: int) -> Iterator[TabulatedBase]:
        """Yield the rows of an indexed table on `refname` from 0-based `start` to `end`.

        The table is BGZF with a tabix or CSI index beside it, as `streampile tabulate --index`
        writes. The region is half-open, as in BED, so it holds the rows with `pos` from
        `start + 1` to `end`.
        """
        decoder = cls(StringIO(), header=False)
        with pybgzf.IndexedReader(path) as index:
            yield from map(decoder.decode, index.query(refname, start, end))


class TabulationWriter(TsvWriter[TabulatedBase], FixedRecordType):
    """A writer of tabulated bases, to a TSV that describes itself.

    A table opens with `##key=value` metadata lines: its format version, under
    `streampile-tabulation`, the streampile version that wrote it, under `streampile-version`,
    and any metadata given, such as the parameters it was tabulated with. A header line follows,
    whose first column is `#contig`, then one row per base. The metadata and the header are
    written once, before the first row or comment, or when `write_header` is first called.

    Attributes:
        metadata: the metadata lines to write, in order.
    """

    @override
    def __init__(
        self,
        handle: TextIO,
        /,
        *,
        metadata: Mapping[str, object] | None = None,
        **options: Unpack[WriterOptions],
    ) -> None:
        """Start a writer with the tabulation's codecs and header unless others are given.

        Args:
            handle: a file-like object to write the table to.
            metadata: more `##key=value` lines to write after the versions, in order.
            options: the options of the writer, with tabulation defaults for any not given.
        """
        _ = options.setdefault("codecs", TABULATION_CODECS)
        _ = options.setdefault("quoting", False)
        _ = options.setdefault("rename", TABULATION_RENAME)
        _ = options.setdefault("comment_prefixes", (METADATA_PREFIX,))
        super().__init__(handle, **options)
        self.metadata: dict[str, str] = {
            "streampile-tabulation": TABULATION_FORMAT,
            "streampile-version": STREAMPILE_VERSION,
            **{key: str(value) for key, value in (metadata or {}).items()},
        }
        self._headed: bool = False

    @override
    def write_header(self) -> None:
        """Write the metadata lines and the header line, unless they are written already."""
        if self._headed:
            return
        self._headed = True
        for key, value in self.metadata.items():
            super().write_comment(f"{METADATA_PREFIX}{key}={value}")
        super().write_header()

    @override
    def write(self, record: TabulatedBase) -> None:
        """Write a row, after the metadata and header if they are not written yet."""
        self.write_header()
        super().write(record)

    @override
    def write_comment(self, comment: str | Comment) -> None:
        """Write a comment, after the metadata and header if they are not written yet."""
        self.write_header()
        super().write_comment(comment)

    @SubscriptableClassmethod
    @classmethod
    def from_path(  # pyright: ignore[reportIncompatibleVariableOverride]
        cls,
        path: Path | str,
        /,
        *,
        index: IndexFormat | None = None,
        index_path: Path | str | None = None,
        threads: int = 1,
        metadata: Mapping[str, object] | None = None,
        **options: Unpack[WriterOptions],
    ) -> Self:
        """Construct a writer of tabulated bases from a file path.

        A path ending in `.gz`, `.bgz`, or `.bgzf` is written as BGZF, which any gzip reader can
        read, and can be indexed with tabix or CSI as it is written, on as many threads as given.
        Rows must then be sorted by position within each contig, and each contig must be contiguous,
        as `tabulate` writes them. The index skips the metadata and header lines, which start with
        `#`. Other paths are written as UTF-8, and compressed when they end in `.bz2` or `.xz`. The
        writer is checked before the file is opened, so a refused writer leaves a file alone.

        Args:
            path: the path to the file to write the table to.
            index: the kind of index to write beside a BGZF file, or None to write none.
            index_path: where to write the index, instead of beside the file; required when the
                file is not a regular file, such as a FIFO.
            threads: the number of threads compressing a BGZF file.
            metadata: more `##key=value` lines to write after the versions, in order.
            options: the options of the writer, with tabulation defaults for any not given.
        """
        path = Path(path).expanduser()
        if path.suffix not in BGZF_SUFFIXES:
            if index is not None or index_path is not None or threads != 1:
                raise ValueError(
                    "An index and threads need a BGZF path ending in .gz, .bgz, or .bgzf, not:"
                    + f" {path}"
                )
            plain: Self = unwrap(super().from_path)(cls, path, metadata=metadata, **options)
            return plain
        if index is None and index_path is not None:
            raise ValueError(f"An index_path needs an index, but none was asked for: {index_path}")
        _ = cls(StringIO(), **options)
        columns = TABULATION_COLUMNS if index is not None else None
        handle = pybgzf.writer(
            path, columns=columns, index=index, index_path=index_path, newline="", threads=threads
        )
        try:
            return cls(handle, metadata=metadata, **options)
        except BaseException:
            handle.close()
            raise
