from array import array
from collections.abc import Iterable
from dataclasses import dataclass
from enum import StrEnum
from enum import auto
from typing import Final
from typing import NamedTuple
from typing import Self
from typing import cast

from pysam import AlignedSegment

from streampile._footprint import DELETED_AT_END
from streampile._footprint import SKIPPED
from streampile._footprint import Footprint
from streampile._footprint import is_placed

DEFAULT_MIN_BASE_QUALITY: Final[int] = 13
"""The default minimum base quality of a pileup, the same as pysam's."""


def _qualities(record: AlignedSegment) -> "array[int] | None":
    return cast("array[int] | None", record.query_qualities)


class PileupReadType(StrEnum):
    """What a read holds at a pileup's position."""

    base = auto()
    deletion = auto()
    insertion = auto()


class PileupRead(NamedTuple):
    """One read at one pileup position.

    A read holding a base or a deletion at a position appears once. A read with an insertion
    right after the position appears again as an insertion entry, as does a read whose alignment
    opens with an insertion, at the position before its first aligned base.

    Attributes:
        alignment: the read.
        query_position: the 0-based query offset of the read's base at the position, or `None`
            for a deletion or an insertion.
        query_position_or_next: the query offset of the read's base at the position, or of the
            read's next base for a deletion, as htslib reports it, or `None` for an insertion or
            a deletion no base follows.
        pileup_type: whether the read holds a base, a deletion, or an insertion.
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
        """The base quality at the position, or of the read's next base for a deletion."""
        qualities = _qualities(self.alignment)
        if self.query_position_or_next is None or qualities is None:
            return None
        return qualities[self.query_position_or_next]

    @property
    def is_del(self) -> bool:
        """Whether the read has a deletion at the position."""
        return self.pileup_type is PileupReadType.deletion

    @property
    def is_ins(self) -> bool:
        """Whether this entry is an insertion right after the position."""
        return self.pileup_type is PileupReadType.insertion

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
        """The base qualities of the inserted bases of an insertion entry, or `None`."""
        qualities = _qualities(self.alignment)
        if self.insertion_offset is None or qualities is None:
            return None
        return list(
            qualities[self.insertion_offset : self.insertion_offset + self.insertion_length]
        )


def pileup_reads(footprint: Footprint, pos: int) -> list[PileupRead]:
    """The entries of one read at one position, which it may not cover."""
    record = footprint.record
    entries: list[PileupRead] = []
    index = pos - footprint.start
    if 0 <= index < len(footprint.offsets):
        offset = footprint.offsets[index]
        if offset >= 0:
            entries.append(PileupRead(record, offset, offset, PileupReadType.base))
        elif offset == DELETED_AT_END:
            entries.append(PileupRead(record, None, None, PileupReadType.deletion))
        elif offset != SKIPPED:
            entries.append(PileupRead(record, None, -offset - 3, PileupReadType.deletion))
    if footprint.insertions:
        insertion = footprint.insertions.get(pos)
        if insertion is not None:
            entries.append(
                PileupRead(
                    record,
                    None,
                    None,
                    PileupReadType.insertion,
                    insertion_offset=insertion[0],
                    insertion_length=insertion[1],
                )
            )
    return entries


@dataclass(frozen=True)
class Pileup:
    """The reads at one reference position.

    Reads with an insertion right after the position are included as insertion entries, so one
    read can have two entries.

    Attributes:
        reference_name: the name of the contig.
        reference_pos: the 0-based position on the contig.
        pileups: the entries of the reads at this position.
        min_base_quality: the base quality below which bases are left out of `filtered_depth`,
            `get_query_qualities`, and `get_query_sequences`.
    """

    reference_name: str
    reference_pos: int
    pileups: list[PileupRead]
    min_base_quality: int = DEFAULT_MIN_BASE_QUALITY

    @property
    def unfiltered_depth(self) -> int:
        """The number of reads with a base or a deletion at this position, ignoring quality."""
        return sum(1 for read in self.pileups if read.pileup_type is not PileupReadType.insertion)

    @property
    def filtered_depth(self) -> int:
        """The number of reads with a base or a deletion at this position at the quality floor.

        A deletion is judged by the quality of the read's next base, as pysam does.
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
    def get_query_qualities(self) -> list[int]:
        """The base qualities of the bases at this position at the quality floor."""
        floor = self.min_base_quality
        return [
            quality
            for read in self.pileups
            if read.pileup_type is PileupReadType.base
            and read.base is not None
            and (quality := read.qual) is not None
            and quality >= floor
        ]

    @property
    def get_query_sequences(self) -> list[str]:
        """The upper-cased bases at this position at the quality floor."""
        floor = self.min_base_quality
        return [
            base
            for read in self.pileups
            if read.pileup_type is PileupReadType.base
            and (base := read.base) is not None
            and (quality := read.qual) is not None
            and quality >= floor
        ]

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
        entries = [
            entry
            for alignment in alignments
            if is_placed(alignment) and alignment.reference_name == contig
            for entry in pileup_reads(Footprint(alignment), pos)
        ]
        return cls(
            reference_name=contig,
            reference_pos=pos,
            pileups=entries,
            min_base_quality=min_base_quality,
        )
