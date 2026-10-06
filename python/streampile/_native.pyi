from collections.abc import Callable
from collections.abc import Iterable
from collections.abc import Iterator
from types import TracebackType
from typing import Any
from typing import ClassVar
from typing import Final
from typing import final

from pysam import AlignedSegment
from pysam import AlignmentFile
from pysam import AlignmentHeader
from pysam import FastaFile
from typing_extensions import Self
from typing_extensions import override

from streampile._pileup import PileupReadType
from streampile._table import TabulatedBase

DEFAULT_EXCLUDE_FLAGS: Final[int]
DEFAULT_MIN_BASE_QUALITY: Final[int]

def direct_bridge(enabled: bool | None = None) -> bool:
    """Whether records are read from htslib's `bam1_t`, and with `enabled`, set it first."""

def normalize_allele(
    pos: int,
    reference_bases: str,
    alternate_bases: str,
    reference: Callable[[int, int], str],
    floor: int = 0,
) -> tuple[int, str, str] | None:
    """Trim and left-align an allele as `bcftools norm` does."""

@final
class Tabulation:
    """The per-read and per-base work of a `Tabulator`, over one reference."""

    def __new__(
        cls,
        reference: FastaFile,
        *,
        min_base_quality: int = 0,
        min_mapping_quality: int = 0,
        exclude_flags: int = ...,
    ) -> Self: ...
    def accepts(self, record: AlignedSegment) -> bool:
        """Whether a read passes the flag and mapping-quality filters and has bases to count."""
    def alleles(
        self, record: AlignedSegment
    ) -> tuple[list[tuple[int, str, str, int, int]], list[tuple[int, int]]]:
        """The alleles a read is counted for, and the stretches of those it is not counted for."""
    def tabulate_contig(
        self, alignments: AlignmentFile, contig: str, spans: list[tuple[int, int]]
    ) -> Iterator[TabulatedBase]:
        """The rows of every base of the spans of one contig, in order."""

@final
class PileupRead:
    """One read at one pileup position."""

    _fields: ClassVar[tuple[str, str, str, str, str, str]]

    def __new__(
        cls,
        alignment: AlignedSegment,
        query_position: int | None,
        query_position_or_next: int | None,
        pileup_type: PileupReadType,
        insertion_offset: int | None = None,
        insertion_length: int = 0,
    ) -> Self: ...
    @property
    def alignment(self) -> AlignedSegment:
        """The read."""
    @property
    def query_position(self) -> int | None:
        """The 0-based query offset of the read's base at the position, or `None` without one."""
    @property
    def query_position_or_next(self) -> int | None:
        """The query offset of the read's base here, or of its next base for a deletion."""
    @property
    def pileup_type(self) -> PileupReadType:
        """Whether the read holds a base, a deletion, a skip, or an insertion."""
    @property
    def insertion_offset(self) -> int | None:
        """The query offset of the first inserted base of an insertion entry."""
    @property
    def insertion_length(self) -> int:
        """The number of inserted bases of an insertion entry."""
    @property
    def base(self) -> str | None:
        """The upper-cased read base at the position, or `None` without one."""
    @property
    def qual(self) -> int | None:
        """The base quality at the position, or of the read's next base for a deletion."""
    @property
    def is_del(self) -> bool:
        """Whether the read has a deletion at the position."""
    @property
    def is_no_call(self) -> bool:
        """Whether the read holds an `N` base, a no-call, at the position."""
    @property
    def is_ins(self) -> bool:
        """Whether this entry is an insertion right after the position."""
    @property
    def is_refskip(self) -> bool:
        """Whether the read skips over the position, with the CIGAR `N` operator."""
    @property
    def inserted_bases(self) -> str | None:
        """The upper-cased inserted bases of an insertion entry, or `None` for any other."""
    @property
    def inserted_qualities(self) -> list[int] | None:
        """The base qualities of the inserted bases of an insertion entry, or `None`."""
    @property
    def five_prime_distance(self) -> int | None:
        """The distance of the read's base from its 5′ end, in bases as sequenced."""
    @property
    def template_end_distance(self) -> int | None:
        """The template's bases between the position and the 5′ end of the read's mate in an FR pair."""
    def _asdict(self) -> dict[str, Any]: ...
    def _replace(
        self,
        *,
        alignment: AlignedSegment = ...,
        query_position: int | None = ...,
        query_position_or_next: int | None = ...,
        pileup_type: PileupReadType = ...,
        insertion_offset: int | None = ...,
        insertion_length: int = ...,
    ) -> Self: ...
    def __len__(self) -> int: ...
    def __getitem__(self, index: int) -> Any: ...
    def __iter__(self) -> Iterator[Any]: ...
    @override
    def __hash__(self) -> int: ...
    @override
    def __eq__(self, other: object) -> bool: ...
    def __lt__(self, other: PileupRead | tuple[Any, ...]) -> bool: ...
    def __le__(self, other: PileupRead | tuple[Any, ...]) -> bool: ...
    def __gt__(self, other: PileupRead | tuple[Any, ...]) -> bool: ...
    def __ge__(self, other: PileupRead | tuple[Any, ...]) -> bool: ...

@final
class Pileup:
    """The reads at one reference position."""

    def __new__(
        cls,
        reference_name: str,
        reference_pos: int,
        pileups: Iterable[PileupRead],
        min_base_quality: int = ...,
    ) -> Self: ...
    @classmethod
    def from_alignments(
        cls,
        alignments: Iterable[AlignedSegment],
        contig: str,
        pos: int,
        min_base_quality: int = ...,
    ) -> Self:
        """Build the pileup at one position from reads in any order."""
    @property
    def reference_name(self) -> str:
        """The name of the contig."""
    @property
    def reference_pos(self) -> int:
        """The 0-based position on the contig."""
    @property
    def pileups(self) -> tuple[PileupRead, ...]:
        """The entries of the reads at this position."""
    @property
    def min_base_quality(self) -> int:
        """The base quality below which bases are left out of the filtered views."""
    @property
    def unfiltered_depth(self) -> int:
        """The number of reads with a base, a deletion, or a skip at this position."""
    @property
    def filtered_depth(self) -> int:
        """The number of reads with a base or a deletion at this position at the quality floor."""
    @property
    def bases(self) -> list[str]:
        """The upper-cased bases at this position at the quality floor, in the order of `pileups`."""
    @property
    def qualities(self) -> list[int]:
        """The base qualities of the bases at this position at the quality floor."""
    def without_overlaps(self) -> Self:
        """A copy of this pileup with one read per template, by query name."""
    @override
    def __hash__(self) -> int: ...
    @override
    def __eq__(self, other: object) -> bool: ...

@final
class StreamingPileupBuilder:
    """Build pileups from coordinate-sorted reads in one forward pass."""

    def __new__(
        cls,
        records: Iterable[AlignedSegment],
        *,
        min_mapping_quality: int = 0,
        exclude_flags: int = ...,
        min_base_quality: int = ...,
        proper_pairs_only: bool = False,
        read_filter: Callable[[AlignedSegment], bool] | None = None,
        tap: Callable[[AlignedSegment], Any] | None = None,
    ) -> Self: ...
    def __enter__(self) -> Self: ...
    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None: ...
    def close(self) -> None:
        """Stop, first handing every read not yet handed to `tap` to it, in input order."""
    def accepts(self, record: AlignedSegment) -> bool:
        """Whether a read passes the built-in filters, is placed, and then passes `read_filter`."""
    def pileup(self, contig: str, pos: int) -> Pileup:
        """Advance to a position, at or after the last one, and pile up the reads there."""
    def columns(self, contig: str, start: int, end: int) -> Iterator[Pileup]:
        """Yield the pileup at every position from `start` to `end`, covered or not."""
    @property
    def header(self) -> AlignmentHeader | None:
        """The header of the reads, or `None` without one."""
    @property
    def min_mapping_quality(self) -> int:
        """The lowest mapping quality of a read to pile up."""
    @property
    def exclude_flags(self) -> int:
        """The SAM flags of reads that are not piled up."""
    @property
    def min_base_quality(self) -> int:
        """The quality floor of each pileup's filtered views."""
    @property
    def proper_pairs_only(self) -> bool:
        """Whether only reads flagged as in a proper pair are piled up."""
    @property
    def read_filter(self) -> Callable[[AlignedSegment], bool] | None:
        """The function that keeps a read for pileups, if any."""
    @property
    def tap(self) -> Callable[[AlignedSegment], Any] | None:
        """The function given every read once the builder has moved past it, if any."""
    @property
    def previous_pileup(self) -> Pileup | None:
        """The last pileup built, which a repeated position returns again."""
