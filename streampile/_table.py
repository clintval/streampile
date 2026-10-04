from dataclasses import dataclass
from inspect import unwrap
from io import StringIO
from pathlib import Path
from typing import Final
from typing import TextIO

import pybgzf
from pybgzf import Columns
from pybgzf import IndexFormat
from typeline import Codecs
from typeline import Comment
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
    spanning the base. A read that skips over the base with an `N` operator observed no base
    there, so it is neither a reference read nor a read of an allele, and is not in `depth`,
    unlike in a pileup's `unfiltered_depth`.

    The alleles anchored at the base are normalized VCF alleles at `pos`, held in parallel
    tuples, most reads first: allele `i` has reference bases `alt_refs[i]`, alternate bases
    `alts[i]`, and `alt_reads[i]` reads, e.g. `C`, `T` for an SNV, `CA`, `C` for a deletion, or
    `C`, `CT` for an insertion.

    Attributes:
        contig: the name of the contig.
        pos: the 1-based position of the base, as in VCF.
        ref: the upper-cased reference base.
        depth: the number of informative reads.
        ref_reads: the number of informative reads with no allele at or spanning the base.
        alt_refs: the reference bases of each allele.
        alts: the alternate bases of each allele.
        alt_reads: the number of reads of each allele.
    """

    contig: str
    pos: int
    ref: str
    depth: int
    ref_reads: int
    alt_refs: tuple[str, ...] = ()
    alts: tuple[str, ...] = ()
    alt_reads: tuple[int, ...] = ()


TABULATION_CODECS: Final[Codecs] = {
    tuple[str, ...]: delimited(str, container=tuple),
    tuple[int, ...]: delimited(int, container=tuple),
}
"""How the allele fields of a `TabulatedBase` are read and written, comma-separated."""

TABULATION_COLUMNS: Final[Columns] = Columns(
    refname=1, start=2, end=None, zero_based=False, meta_char="#", skip_lines=1
)
"""Where an index finds each row: the contig and 1-based position, after the header line."""

BGZF_SUFFIXES: Final[tuple[str, ...]] = (".bgz", ".gz")
"""The file extensions written as BGZF, which is also valid gzip."""


class TabulationReader(TsvReader[TabulatedBase], FixedRecordType):
    """A reader of tabulated bases, from a TSV with a header, gzipped or BGZF or not."""

    @override
    def __init__(self, handle: TextIO, /, **options: Unpack[ReaderOptions]) -> None:
        """Start a reader with the tabulation's codecs unless others are given.

        Args:
            handle: a file-like object to read the table from.
            options: the options of the reader, with tabulation defaults for any not given.
        """
        _ = options.setdefault("codecs", TABULATION_CODECS)
        _ = options.setdefault("quoting", False)
        super().__init__(handle, **options)


class TabulationWriter(TsvWriter[TabulatedBase], FixedRecordType):
    """A writer of tabulated bases, to a TSV."""

    @override
    def __init__(self, handle: TextIO, /, **options: Unpack[WriterOptions]) -> None:
        """Start a writer with the tabulation's codecs unless others are given.

        Args:
            handle: a file-like object to write the table to.
            options: the options of the writer, with tabulation defaults for any not given.
        """
        _ = options.setdefault("codecs", TABULATION_CODECS)
        _ = options.setdefault("quoting", False)
        super().__init__(handle, **options)
        self._indexed: bool = False
        self._headed: bool = False

    @override
    def write_header(self) -> None:
        """Write the header line."""
        super().write_header()
        self._headed = True

    @override
    def write(self, record: TabulatedBase) -> None:
        """Write a row, refusing it before the header when the table is being indexed."""
        self._check_headed("a row")
        super().write(record)

    @override
    def write_comment(self, comment: str | Comment) -> None:
        """Write a comment, refusing it before the header when the table is being indexed."""
        self._check_headed("a comment")
        super().write_comment(comment)

    def _check_headed(self, what: str) -> None:
        """Refuse a line before the header of an indexed table."""
        if self._indexed and not self._headed:
            raise ValueError(
                f"Cannot write {what} to an indexed table before its header, which the index"
                + " skips as the first line! Call write_header() first."
            )

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
        **options: Unpack[WriterOptions],
    ) -> Self:
        """Construct a writer of tabulated bases from a file path.

        A path ending in `.gz` or `.bgz` is written as BGZF, which any gzip reader can read, and
        can be indexed with tabix or CSI as it is written, on as many threads as given.
        Rows must then be sorted by position within each contig, and each contig must be
        contiguous, as `tabulate` writes them. The header must also come first, since the index
        skips the first line: without it the index would skip the first row, and after a comment
        it would read the header as a row, so rows and comments are refused until it is written.
        Other paths are written as UTF-8, and compressed when they end in `.bz2` or `.xz`.
        The writer is checked before the file is opened, so a refused writer leaves a file alone.

        Args:
            path: the path to the file to write the table to.
            index: the kind of index to write beside a BGZF file, or None to write none.
            index_path: where to write the index, instead of beside the file; required when the
                file is not a regular file, such as a FIFO.
            threads: the number of threads compressing a BGZF file.
            options: the options of the writer, with tabulation defaults for any not given.
        """
        path = Path(path).expanduser()
        if path.suffix not in BGZF_SUFFIXES:
            if index is not None or index_path is not None or threads != 1:
                raise ValueError(
                    f"An index and threads need a BGZF path ending in .gz or .bgz, not: {path}"
                )
            plain: Self = unwrap(super().from_path)(cls, path, **options)
            return plain
        if index is None and index_path is not None:
            raise ValueError(f"An index_path needs an index, but none was asked for: {index_path}")
        _ = cls(StringIO(), **options)
        columns = TABULATION_COLUMNS if index is not None else None
        handle = pybgzf.writer(
            path, columns=columns, index=index, index_path=index_path, newline="", threads=threads
        )
        try:
            writer = cls(handle, **options)
        except BaseException:
            handle.close()
            raise
        writer._indexed = index is not None
        return writer
