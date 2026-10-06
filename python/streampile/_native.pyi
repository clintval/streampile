from collections.abc import Callable
from collections.abc import Iterable
from collections.abc import Iterator
from dataclasses import Field
from types import TracebackType
from typing import Any
from typing import ClassVar
from typing import Final
from typing import SupportsIndex
from typing import final

from pysam import AlignedSegment
from pysam import AlignmentFile
from pysam import AlignmentHeader
from pysam import FastaFile
from typing_extensions import Self
from typing_extensions import override

from streampile._pileup import AgreementStrategy
from streampile._pileup import DisagreementStrategy
from streampile._pileup import PileupReadType
from streampile._table import TabulatedBase

DEFAULT_EXCLUDE_FLAGS: Final[int]
DEFAULT_MIN_BASE_QUALITY: Final[int]

def direct_bridge(enabled: bool | None = None) -> bool:
    """Whether records are read from htslib's `bam1_t`, and with `enabled`, sets it first: records
    are read from `bam1_t` only when its layout was found as expected.
    """

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
    @property
    def reference(self) -> FastaFile:
        """The indexed reference the reads are aligned to."""
    @property
    def min_base_quality(self) -> int:
        """The lowest base quality of an informative base."""
    @property
    def min_mapping_quality(self) -> int:
        """The lowest mapping quality of a counted read."""
    @property
    def exclude_flags(self) -> int:
        """The SAM flags of reads that are not counted."""
    def accepts(self, record: AlignedSegment) -> bool:
        """Whether a read passes the flag and mapping-quality filters and has bases to count."""
    def alleles(
        self, record: AlignedSegment
    ) -> tuple[list[tuple[int, str, str, int, int]], list[tuple[int, int]]]:
        """The alleles a read is counted for, as `(pos, ref, alt, start, end)`, and the stretches of
        those it is not counted for.
        """
    def tabulate_contig(
        self, alignments: AlignmentFile, contig: str, spans: list[tuple[int, int]]
    ) -> ContigRows:
        """The rows of every base of the spans of one contig, in order, from reads fetched span by
        span from an indexed alignment file.
        """

@final
class ContigRows:
    """The rows of every base of the spans of one contig, made one at a time."""

    def __iter__(self) -> Self: ...
    def __next__(self) -> TabulatedBase: ...

@final
class PileupRead:
    """One read at one pileup position.

    A read holding a base, a deletion, or a reference skip (the CIGAR `N` operator) at a position
    appears once. A read with an insertion right after the position appears again as an insertion
    entry, as does a read whose alignment opens with an insertion, at the position before its first
    aligned base. So an insertion at either end of an alignment is reported: htslib reports only one
    that closes an alignment, and fgbio only one that opens it, at offset 0 even after a soft clip.

    A skip entry holds no base, no quality, and no query offset. htslib flags the same entry as
    both `is_del` and `is_refskip` and gives it the offset and quality of the read's next base;
    here a skip is not a deletion, and has no quality to pass a floor with.

    A `PileupRead` behaves as the `NamedTuple` of its six fields, though it is not a `tuple`: it
    unpacks, indexes, compares, hashes, matches, copies, and pickles as one, and has `_make`,
    `_asdict`, and `_replace`. Its bases and qualities are those of the read as it was when it was
    piled up, while `alignment` is the very object given to the builder, so it can be changed,
    e.g. tagged, and written on.

    Attributes:
        alignment: the read.
        query_position: the 0-based query offset of the read's base at the position, or `None`
            for a deletion, a skip, or an insertion.
        query_position_or_next: the query offset of the read's base at the position, or of the
            read's next base for a deletion, as htslib reports it, or `None` for a skip, an
            insertion, or a deletion no base follows.
        pileup_type: whether the read holds a base, a deletion, a skip, or an insertion.
        insertion_offset: the query offset of the first inserted base of an insertion entry.
        insertion_length: the number of inserted bases of an insertion entry.
    """

    _fields: ClassVar[tuple[str, str, str, str, str, str]]
    _field_defaults: ClassVar[dict[str, Any]]
    __match_args__ = (
        "alignment",
        "query_position",
        "query_position_or_next",
        "pileup_type",
        "insertion_offset",
        "insertion_length",
    )

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
        """The query offset of the read's base here, or of its next base for a deletion, or `None`."""
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
        """The upper-cased read base at the position, or `None` without one.

        A base written as `=`, a match to the reference, is `=`, as a pileup knows no reference.
        """
    @property
    def qual(self) -> int | None:
        """The base quality at the position, or of the read's next base for a deletion.

        A read with bases but no stored qualities (QUAL `*`) has quality 255 at every base, as in
        htslib, so it passes every floor. It is `None` where there is no base to take it from: for
        an insertion, a deletion no base follows, or a read with no stored bases (SEQ `*`).
        """
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
        """The base qualities of the inserted bases of an insertion entry, or `None`.

        They are 255 for a read with no stored qualities, and `None` for one with no stored bases.
        """
    @property
    def five_prime_distance(self) -> int | None:
        """The distance of the read's base from its 5′ end, in bases as sequenced.

        It is the query offset for a forward read, counted from the other end for a reverse read,
        so soft-clipped bases count, and 0 is the first base sequenced: fgbio's
        `positionInReadInReadOrder` minus one. For a deletion or a skip, which holds no base, it is
        the number of the read's bases sequenced before the position, as `template_end_distance`
        counts the bases after it. It is `None` for an insertion entry, and for an entry made by hand
        that holds no base, whose position is unknown.
        """
    @property
    def template_end_distance(self) -> int | None:
        """The number of the template's bases between the position and the template's other end, the
        5′ end of the mate of a read in an FR pair: 0 at the mate's 5′ end.

        It walks the read's CIGAR and the mate's, from its `MC` tag, so an indel counts by its
        length; soft-clipped bases count and hard-clipped bases, absent from the records, do not.
        Where the mate has no base, a position both reads align to carries the count, or else the
        reference between them; the template length (TLEN) is never read. It is `None` for a read
        that `is_fr_pair` does not call a read of an FR pair, for a read of an FR pair only at a
        position past the mate's 5′ end, where a read runs through its mate, and for an entry made
        by hand that holds no base, whose position is unknown.

        Raises:
            ValueError: for a forward read whose mate is mapped to its contig on the other strand,
                or a reverse read of an FR pair, with no `MC` tag, or one that is not a CIGAR
                string, or, for a read of an FR pair, one that spans no reference.
        """
    @property
    def is_fr_pair(self) -> bool:
        """Whether the read is a read of an FR pair, as htsjdk 5.0.0's `getPairOrientation` says.

        A pair is FR when its reads are paired, mapped to one contig, and on opposite strands, and
        the forward read's aligned 5′ position is at or before the reverse read's, so a pair whose
        5′ ends coincide is FR. A forward read takes its mate's aligned end from its `MC` tag and
        never from the template length (TLEN).

        Raises:
            ValueError: for a forward read whose mate is mapped to its contig on the other strand,
                whatever the pair's orientation, with no `MC` tag, or one that is not a CIGAR string.
        """
    @classmethod
    def _make(cls, iterable: Iterable[Any]) -> Self: ...
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
    def __replace__(
        self,
        *,
        alignment: AlignedSegment = ...,
        query_position: int | None = ...,
        query_position_or_next: int | None = ...,
        pileup_type: PileupReadType = ...,
        insertion_offset: int | None = ...,
        insertion_length: int = ...,
    ) -> Self: ...
    def count(self, value: Any, /) -> int: ...
    def index(
        self, value: Any, start: SupportsIndex = 0, stop: SupportsIndex = ..., /
    ) -> int: ...
    def __contains__(self, value: object, /) -> bool: ...
    def __add__(self, other: tuple[Any, ...] | PileupRead, /) -> tuple[Any, ...]: ...
    def __radd__(self, other: tuple[Any, ...], /) -> tuple[Any, ...]: ...
    def __mul__(self, times: SupportsIndex, /) -> tuple[Any, ...]: ...
    def __rmul__(self, times: SupportsIndex, /) -> tuple[Any, ...]: ...
    def __getnewargs__(self) -> tuple[Any, ...]: ...
    def __copy__(self) -> Self: ...
    def __deepcopy__(self, memo: dict[int, Any], /) -> Self: ...
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
class PileupTemplate:
    """One template at one pileup position: the reads of one query name, their bases called into one.

    A template's strand and distances are those of its first read, the first of a pair or a
    fragment's only read, worked out from its second read where the first holds no base here.

    Attributes:
        query_name: the name of the template's reads, `*` for a read with none, which is a
            template of its own.
        reads: the entries of the template's reads at the position, usually one or two.
        pileup_type: whether the template holds a base, a deletion, or a skip.
        base: the template's base, its reads' bases called into one.
        qual: the quality of the template's base, or of the next base for a deletion.
    """

    @property
    def query_name(self) -> str:
        """The name of the template's reads, `*` for a read with none, which is a template of its own."""
    @property
    def reads(self) -> tuple[PileupRead, ...]:
        """The entries of the template's reads at the position, usually one or two, as in `pileups`."""
    @property
    def pileup_type(self) -> PileupReadType:
        """Whether the template holds a base, a deletion, or a skip: a base if a read's base votes, or
        else a deletion if a read's deletion votes, or else, with no vote, what its reads hold.
        """
    @property
    def base(self) -> str | None:
        """The template's upper-cased base, its voting reads' bases called into one, or `None`."""
    @property
    def qual(self) -> int | None:
        """The quality of the template's base, or of the next base for a deletion, or `None` without
        a vote.
        """
    @property
    def is_del(self) -> bool:
        """Whether the template holds a deletion at the position, and no base that votes."""
    @property
    def is_refskip(self) -> bool:
        """Whether every read of the template skips over the position."""
    @property
    def is_no_call(self) -> bool:
        """Whether the template's base is a no-call, `N`."""
    @property
    def is_reverse(self) -> bool:
        """Whether the template's first read is aligned to the reverse strand.

        It is `False` for an F1R2 pair and `True` for an F2R1 pair, read from the second read's
        flags without the first.
        """
    @property
    def five_prime_distance(self) -> int | None:
        """The number of the template's bases between its first read's 5′ end and the position.

        It is the first read's `five_prime_distance` where it holds a base here, and otherwise the
        second read's `template_end_distance`.

        Raises:
            ValueError: where the read's `template_end_distance` or `is_fr_pair` raises.
        """
    @property
    def template_end_distance(self) -> int | None:
        """The number of the template's bases between the position and its other end, the 5′ end of
        the second read of an FR pair.

        It is the first read's `template_end_distance`, and without the first read here the second
        read's `five_prime_distance` where that read `is_fr_pair`; it is `None` for any pair that is
        not FR.

        Raises:
            ValueError: where the read's `template_end_distance` or `is_fr_pair` raises.
        """

@final
class Pileup:
    """The reads at one reference position.

    Reads with an insertion right after the position are included as insertion entries, so one read
    can have two entries. Reads that skip over the position with the CIGAR `N` operator are included
    as skip entries. A pileup is a snapshot: it outlives the builder's next move, and its views are
    worked out from it in Rust. It behaves as a frozen dataclass of its four fields: it compares,
    hashes, matches, copies, and pickles by them, and `dataclasses.replace`, `fields`, `asdict`,
    and `astuple` take it.

    Attributes:
        reference_name: the name of the contig.
        reference_pos: the 0-based position on the contig.
        pileups: the entries of the reads at this position.
        min_base_quality: the base quality below which bases are left out of `filtered_depth`,
            `bases`, `qualities`, and the votes of `templates()`.
    """

    __match_args__ = ("reference_name", "reference_pos", "pileups", "min_base_quality")
    __dataclass_fields__: ClassVar[dict[str, Field[Any]]]

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
        """Build the pileup at one position from reads in any order.

        Unmapped reads and reads on other contigs are ignored.

        Args:
            alignments: the reads to pile up.
            contig: the name of the contig.
            pos: the 0-based position on the contig.
            min_base_quality: the quality floor of the pileup's filtered views and template votes.
        """
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
        """The base quality below which bases are left out of the filtered views and template votes."""
    @property
    def unfiltered_depth(self) -> int:
        """The number of reads with a base, a deletion, or a skip at this position.

        Quality is ignored, and insertion entries are not counted, as in htslib's column depth.
        """
    @property
    def filtered_depth(self) -> int:
        """The number of reads with a base or a deletion at this position at the quality floor.

        A deletion is judged by the quality of the read's next base, as pysam does. A skip has no
        quality, so it is never counted, while pysam judges a skip by its next base too.
        """
    @property
    def bases(self) -> list[str]:
        """The upper-cased bases at this position at the quality floor, in the order of `pileups`.

        Only base entries at `min_base_quality` are listed, so the list lines up with `qualities`.
        pysam's `get_query_sequences()`, by contrast, lists every entry, with an empty string for a
        deletion or a skip.
        """
    @property
    def qualities(self) -> list[int]:
        """The base qualities of the bases at this position at the quality floor, as in `bases`.

        pysam's `get_query_qualities()`, by contrast, lists every entry, with the quality of the
        read's next base for a deletion or a skip.
        """
    def templates(
        self,
        *,
        agreement: AgreementStrategy = ...,
        disagreement: DisagreementStrategy = ...,
    ) -> list[PileupTemplate]:
        """One observation per template at this position, its reads grouped by query name, in the
        order of each template's first entry in `pileups`. A read with no name is a template of its
        own.

        Where two reads of a template hold bases, they are called into one as fgbio's
        `CallOverlappingConsensusBases` calls them, as fgumi implements it: `agreement` makes the
        quality of equal bases and `disagreement` the base and quality of different ones. Where the
        strategies leave a read's base or quality unchanged, the template takes the higher quality.

        fgumi defines no more than that, so at each position a no-call (`N`) is left alone, as fgumi
        leaves it: the other read's base stands at its own quality, and two no-calls are an `N` at
        the higher quality. A no-call is still a base, so it stands over the other read's deletion,
        which fgumi leaves alone too. A read with a deletion or a skip holds no base, so a template
        whose other read holds one has that base at its own quality; a template with no base is a
        deletion if either read holds one, at the higher of their qualities, or else a skip.
        Insertion entries are no part of a template, so a read whose only entry here is an insertion
        adds nothing to its template. A read under `min_base_quality` does not vote, as `bases`,
        `qualities`, and `filtered_depth` leave it out: a base, or a deletion judged by its next
        base, under the floor. So a mate under the floor neither masks nor lowers the other mate's
        base, and a template none of whose reads votes has no `base` or `qual`. A template with more
        than two reads here, as when supplementary records are piled up, calls them in the order of
        `pileups`.

        Args:
            agreement: how the quality of two reads holding the same base is made.
            disagreement: how the base and quality of two reads holding different bases are made.
        """
    def __replace__(
        self,
        *,
        reference_name: str = ...,
        reference_pos: int = ...,
        pileups: Iterable[PileupRead] = ...,
        min_base_quality: int = ...,
    ) -> Self: ...
    def __copy__(self) -> Self: ...
    def __deepcopy__(self, memo: dict[int, Any], /) -> Self: ...
    @override
    def __hash__(self) -> int: ...
    @override
    def __eq__(self, other: object) -> bool: ...

@final
class StreamingPileupBuilder:
    """Build pileups from coordinate-sorted reads in one forward pass.

    Ask for pileups at positions that never move backwards: the same position again returns the
    pileup already built, and an earlier one raises a `ValueError`. Each read's CIGAR is walked
    once, when the read is first reached, so building a pileup costs one lookup per read there.
    Reads are piled up in Rust, from a copy of each read's record made as it is read.

    Every read, filtered or not, is handed to `tap` exactly once and in input order, as soon as
    the builder has moved past it and every read before it, or when the builder closes. A read
    can therefore be changed, e.g. tagged, while it is in a pileup and then written by `tap`.
    Pileups see each read as it was when the builder read it; changes made after that,
    including in `read_filter`, reach `tap` and `alignment` but not later pileups.
    Keeping input order means a read is held until every read before it has been passed, so a
    long read, e.g. one with a long reference skip, holds back every read that starts within it:
    with a `tap`, buffering behind the longest active read is inherent. Without a `tap`, a read
    is dropped as soon as the builder has moved past it.

    An exception while advancing to a position, from `records`, a malformed read, `read_filter`,
    or `tap`, stops the builder, so a pileup never misses a read: every later pileup raises a
    `ValueError`, and closing still hands every read not yet handed over to `tap`.

    A builder can be used from any thread, one call at a time: a call made while another runs,
    such as from `read_filter` or `tap`, raises a `RuntimeError`.

    ```python
    with (
        AlignmentFile("in.bam", threads=4) as source,
        AlignmentFile("out.bam", "wb", template=source) as sink,
        StreamingPileupBuilder(source, tap=sink.write) as builder,
    ):
        pileup = builder.pileup("chr1", 100)
    ```
    """

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
        """Stop, first handing every read not yet handed to `tap` to it, in input order.

        With a `tap`, the rest of the input is read to the end, so an output written by `tap` is
        complete. Without one, no more of the input is read. A read goes to `tap` once even when
        `tap` raises, and closing again hands over the reads after it.
        """
    def accepts(self, record: AlignedSegment) -> bool:
        """Whether a read passes the built-in filters, is placed, and then passes `read_filter`."""
    def pileup(self, contig: str, pos: int) -> Pileup:
        """Advance to a position, at or after the last one, and pile up the reads there.

        Args:
            contig: the name of the contig.
            pos: the 0-based position on the contig.

        Raises:
            ValueError: if the builder is closed or stopped by an earlier exception, the contig is
                not in the header, the position is negative, or the position is before the last
                one asked for.
        """
    def columns(self, contig: str, start: int, end: int) -> Columns:
        """Yield the pileup at every position from `start` to `end`, covered or not.

        Args:
            contig: the name of the contig.
            start: the 0-based first position.
            end: the 0-based position after the last one.

        Raises:
            ValueError: if `end` is before `start`, or as `pileup` does.
        """
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
        """The quality floor of each pileup's filtered views and template votes."""
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

@final
class Columns:
    """The pileup at every position of a span, from `StreamingPileupBuilder.columns`."""

    def __iter__(self) -> Self: ...
    def __next__(self) -> Pileup: ...
