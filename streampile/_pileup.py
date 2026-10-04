from collections.abc import Iterable
from dataclasses import dataclass
from dataclasses import replace
from enum import StrEnum
from enum import auto
from typing import Final
from typing import NamedTuple
from typing import Self

from pysam import AlignedSegment

from streampile._footprint import DELETED_AT_END
from streampile._footprint import SKIPPED
from streampile._footprint import Footprint
from streampile._footprint import is_placed
from streampile._footprint import query_qualities

DEFAULT_MIN_BASE_QUALITY: Final[int] = 13
"""The default minimum base quality of a pileup, the same as pysam's."""

DEFAULT_EXCLUDE_FLAGS: Final[int] = 0xF00
"""Secondary, QC-fail, duplicate, and supplementary reads, which are left out by default."""

MISSING_BASE_QUALITY: Final[int] = 255
"""The base quality of every base of a read with no stored qualities (QUAL `*`), as in htslib."""


class PileupReadType(StrEnum):
    """What a read holds at a pileup's position."""

    base = auto()
    deletion = auto()
    insertion = auto()
    skip = auto()


BASE = PileupReadType.base
DELETION = PileupReadType.deletion
INSERTION = PileupReadType.insertion
SKIP = PileupReadType.skip


class PileupRead(NamedTuple):
    """One read at one pileup position.

    A read holding a base, a deletion, or a reference skip (an `N` operator) at a position appears
    once. A read with an insertion right after the position appears again as an insertion entry,
    as does a read whose alignment opens with an insertion, at the position before its first
    aligned base. So an insertion at either end of an alignment is reported: htslib reports only
    one that closes an alignment, and fgbio only one that opens it, at offset 0 even after a
    soft clip.

    A skip entry holds no base, no quality, and no query offset. htslib flags the same entry as
    both `is_del` and `is_refskip` and gives it the offset and quality of the read's next base;
    here a skip is not a deletion, and has no quality to pass a floor with.

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

    alignment: AlignedSegment
    query_position: int | None
    query_position_or_next: int | None
    pileup_type: PileupReadType
    insertion_offset: int | None = None
    insertion_length: int = 0

    @property
    def base(self) -> str | None:
        """The upper-cased read base at the position, or `None` without one."""
        sequence = self.alignment.query_sequence
        if self.query_position is None or sequence is None:
            return None
        return sequence[self.query_position].upper()

    @property
    def qual(self) -> int | None:
        """The base quality at the position, or of the read's next base for a deletion.

        A read with bases but no stored qualities (QUAL `*`) has quality 255 at every base, as in
        htslib, so it passes every floor. It is `None` where there is no base to take it from: for
        an insertion, a deletion no base follows, or a read with no stored bases (SEQ `*`).
        """
        offset = self.query_position_or_next
        if offset is None:
            return None
        qualities = query_qualities(self.alignment)
        if qualities is None:
            return MISSING_BASE_QUALITY if offset < self.alignment.query_length else None
        return qualities[offset]

    @property
    def is_del(self) -> bool:
        """Whether the read has a deletion at the position."""
        return self.pileup_type is PileupReadType.deletion

    @property
    def is_ins(self) -> bool:
        """Whether this entry is an insertion right after the position."""
        return self.pileup_type is PileupReadType.insertion

    @property
    def is_refskip(self) -> bool:
        """Whether the read skips over the position with an `N` operator."""
        return self.pileup_type is PileupReadType.skip

    @property
    def inserted_bases(self) -> str | None:
        """The upper-cased inserted bases of an insertion entry, or `None` for any other."""
        sequence = self.alignment.query_sequence
        if self.insertion_offset is None or sequence is None:
            return None
        return sequence[
            self.insertion_offset : self.insertion_offset + self.insertion_length
        ].upper()

    @property
    def inserted_qualities(self) -> list[int] | None:
        """The base qualities of the inserted bases of an insertion entry, or `None`.

        They are 255 for a read with no stored qualities, and `None` for one with no stored bases.
        """
        offset = self.insertion_offset
        if offset is None:
            return None
        qualities = query_qualities(self.alignment)
        if qualities is None:
            if not self.alignment.query_length:
                return None
            return [MISSING_BASE_QUALITY] * self.insertion_length
        return list(qualities[offset : offset + self.insertion_length])


def pileup_entries(footprints: Iterable[Footprint], pos: int) -> tuple[PileupRead, ...]:
    """The entries of reads at one position, in the order of the reads, which may not cover it."""
    entries: list[PileupRead] = []
    append = entries.append
    for footprint in footprints:
        record = footprint.record
        index = pos - footprint.start
        if 0 <= index < len(footprint.offsets):
            offset = footprint.offsets[index]
            if offset >= 0:
                append(PileupRead(record, offset, offset, BASE))
            elif offset < DELETED_AT_END:
                append(PileupRead(record, None, -offset - 3, DELETION))
            elif offset == SKIPPED:
                append(PileupRead(record, None, None, SKIP))
            else:
                append(PileupRead(record, None, None, DELETION))
        if footprint.insertions:
            insertion = footprint.insertions.get(pos)
            if insertion is not None:
                append(PileupRead(record, None, None, INSERTION, insertion[0], insertion[1]))
    return tuple(entries)


@dataclass(frozen=True)
class Pileup:
    """The reads at one reference position.

    Reads with an insertion right after the position are included as insertion entries, so one
    read can have two entries. Reads that skip over the position with an `N` operator are included
    as skip entries.

    Attributes:
        reference_name: the name of the contig.
        reference_pos: the 0-based position on the contig.
        pileups: the entries of the reads at this position.
        min_base_quality: the base quality below which bases are left out of `filtered_depth`,
            `bases`, and `qualities`.
    """

    reference_name: str
    reference_pos: int
    pileups: tuple[PileupRead, ...]
    min_base_quality: int = DEFAULT_MIN_BASE_QUALITY

    @property
    def unfiltered_depth(self) -> int:
        """The number of reads with a base, a deletion, or a skip at this position.

        Quality is ignored, and insertion entries are not counted, as in htslib's column depth.
        """
        return sum(1 for read in self.pileups if read.pileup_type is not PileupReadType.insertion)

    @property
    def filtered_depth(self) -> int:
        """The number of reads with a base or a deletion at this position at the quality floor.

        A deletion is judged by the quality of the read's next base, as pysam does. A skip has no
        quality, so it is never counted, while pysam judges a skip by its next base too.
        """
        floor = self.min_base_quality
        return sum(
            1
            for read in self.pileups
            if read.pileup_type is not PileupReadType.insertion
            and (quality := read.qual) is not None
            and quality >= floor
        )

    @property
    def bases(self) -> list[str]:
        """The upper-cased bases at this position at the quality floor, in the order of `pileups`.

        Only base entries at `min_base_quality` are listed, so the list lines up with `qualities`.
        pysam's `get_query_sequences()`, by contrast, lists every entry, with an empty string for a
        deletion or a skip.
        """
        floor = self.min_base_quality
        return [
            base
            for read in self.pileups
            if read.pileup_type is PileupReadType.base
            and (base := read.base) is not None
            and (quality := read.qual) is not None
            and quality >= floor
        ]

    @property
    def qualities(self) -> list[int]:
        """The base qualities of the bases at this position at the quality floor, as in `bases`.

        pysam's `get_query_qualities()`, by contrast, lists every entry, with the quality of the
        read's next base for a deletion or a skip.
        """
        floor = self.min_base_quality
        return [
            quality
            for read in self.pileups
            if read.pileup_type is PileupReadType.base
            and (quality := read.qual) is not None
            and quality >= floor
        ]

    def without_overlaps(self) -> Self:
        """A copy of this pileup with one read per template, by query name.

        The read kept for a template is the first of its name in `pileups` whose entry here is a
        base or a deletion at `min_base_quality`, or else the first of its name. A builder fills
        `pileups` in input order, so of two passing mates in a coordinate-sorted file, the one
        that starts first is kept. A mate's skip, or its base under the floor, therefore never
        hides the other mate's base, as in htslib, and as in fgbio, whose floor drops a failing
        base before its `withoutOverlaps` keeps the first entry. Every entry of the kept read
        stays, its insertion entry included, where fgbio keeps only the first entry of each
        template.
        """
        floor = self.min_base_quality
        first: dict[str | None, AlignedSegment] = {}
        passing: dict[str | None, AlignedSegment] = {}
        for entry in self.pileups:
            name = entry.alignment.query_name
            _ = first.setdefault(name, entry.alignment)
            if (quality := entry.qual) is not None and quality >= floor:
                _ = passing.setdefault(name, entry.alignment)
        kept = first | passing
        pileups = tuple(
            entry for entry in self.pileups if kept[entry.alignment.query_name] is entry.alignment
        )
        return replace(self, pileups=pileups)

    @classmethod
    def from_alignments(
        cls,
        alignments: Iterable[AlignedSegment],
        contig: str,
        pos: int,
        min_base_quality: int = DEFAULT_MIN_BASE_QUALITY,
    ) -> Self:
        """Build the pileup at one position from reads in any order.

        Unmapped reads and reads on other contigs are ignored.

        Args:
            alignments: the reads to pile up.
            contig: the name of the contig.
            pos: the 0-based position on the contig.
            min_base_quality: the quality floor of the pileup's filtered views.
        """
        footprints = (
            Footprint(alignment)
            for alignment in alignments
            if is_placed(alignment) and alignment.reference_name == contig
        )
        return cls(
            reference_name=contig,
            reference_pos=pos,
            pileups=pileup_entries(footprints, pos),
            min_base_quality=min_base_quality,
        )
