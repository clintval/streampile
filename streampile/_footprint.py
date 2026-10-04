from array import array
from typing import Final
from typing import final

from pysam import CDEL
from pysam import CDIFF
from pysam import CEQUAL
from pysam import CINS
from pysam import CMATCH
from pysam import CREF_SKIP
from pysam import CSOFT_CLIP
from pysam import AlignedSegment

SKIPPED: Final[int] = -1
"""The offset of a reference position that the read skips over with an `N` operator."""

DELETED_AT_END: Final[int] = -2
"""The offset of a deleted reference position that no base of the read follows."""

REFERENCE_OPERATORS: Final[frozenset[int]] = frozenset({CMATCH, CDEL, CREF_SKIP, CEQUAL, CDIFF})
"""The CIGAR operators that consume the reference: M, D, N, =, and X."""


@final
class Footprint:
    """Where one read sits on the reference, worked out once from its CIGAR.

    Attributes:
        record: the read.
        start: the 0-based reference position of the read's first aligned or deleted base.
        end: the 0-based reference position just past the read's last aligned or deleted base.
        offsets: one entry per reference position from `start` to `end`. A non-negative entry is
            the query offset of the base aligned there, `SKIPPED` marks a reference skip, and a
            deleted position holds `DELETED_AT_END` or `-(next) - 3`, where `next` is the query
            offset of the read's next base, as htslib reports it.
        insertions: the query offset and length of each insertion, keyed by the reference
            position just before it.
    """

    __slots__ = ("end", "insertions", "offsets", "record", "start")

    def __init__(self, record: AlignedSegment) -> None:
        """Walk the CIGAR of a mapped read with at least one reference-consuming operator."""
        start: int = record.reference_start
        offsets = array("i")
        insertions: dict[int, tuple[int, int]] = {}
        query_length: int = record.query_length or record.infer_query_length() or 0
        query = 0
        position = start
        for operator, length in record.cigartuples or ():
            if operator == CMATCH or operator == CEQUAL or operator == CDIFF:
                offsets.extend(range(query, query + length))
                query += length
                position += length
            elif operator == CINS:
                anchor = position - 1
                previous = insertions.get(anchor)
                if previous is not None and previous[0] + previous[1] == query:
                    insertions[anchor] = (previous[0], previous[1] + length)
                else:
                    insertions[anchor] = (query, length)
                query += length
            elif operator == CDEL:
                offsets.extend([-query - 3 if query < query_length else DELETED_AT_END] * length)
                position += length
            elif operator == CREF_SKIP:
                offsets.extend([SKIPPED] * length)
                position += length
            elif operator == CSOFT_CLIP:
                query += length
        self.record: AlignedSegment = record
        self.start: int = start
        self.end: int = position
        self.offsets: array[int] = offsets
        self.insertions: dict[int, tuple[int, int]] = insertions


def is_placed(record: AlignedSegment) -> bool:
    """Whether a read is mapped and has a reference-consuming operator: M, D, N, =, or X.

    A read with none, such as `4S`, `4I`, or `2S2I`, is left out, as htslib leaves it out of a
    pileup, although pysam gives it a `reference_end` one past its start.
    """
    return (
        not record.is_unmapped
        and record.reference_id >= 0
        and any(operator in REFERENCE_OPERATORS for operator, _ in record.cigartuples or ())
    )
