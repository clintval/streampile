import re
import sys
from collections import deque
from collections.abc import Callable
from collections.abc import Iterable
from collections.abc import Iterator
from types import TracebackType
from typing import Any
from typing import Self

from pysam import AlignedSegment
from pysam import AlignmentHeader

from streampile._footprint import DELETED_AT_END
from streampile._footprint import SKIPPED
from streampile._footprint import Footprint
from streampile._footprint import is_placed
from streampile._pileup import DEFAULT_MIN_BASE_QUALITY
from streampile._pileup import Pileup
from streampile._pileup import PileupRead
from streampile._pileup import PileupReadType

SORT_ORDER = re.compile(r"^@HD\t.*\bSO:([^\t\n]+)", re.MULTILINE)

UNPLACED: int = sys.maxsize
"""The contig index given to reads with no contig, which sort after every placed read."""

BASE = PileupReadType.base
DELETION = PileupReadType.deletion
INSERTION = PileupReadType.insertion


class StreamingPileupBuilder:
    """Build pileups from coordinate-sorted reads in one forward pass.

    Ask for pileups at positions that never move backwards: the same position again returns the
    pileup already built, and an earlier one raises a `ValueError`. Each read's CIGAR is walked
    once, when the read is first reached, so building a pileup costs one lookup per read there.

    Every read, filtered or not, is handed to `tap` exactly once and in input order, as soon as
    the builder has moved past it and every read before it, or when the builder closes. A read
    can therefore be changed, e.g. tagged, while it is in a pileup and then written by `tap`.
    Keeping input order means a read is held until every read before it has been passed, so a
    long read, e.g. one with a long reference skip, holds back every read that starts within it:
    with a `tap`, buffering behind the longest active read is inherent. Without a `tap`, a read
    is dropped as soon as the builder has moved past it.

    ```python
    with (
        AlignmentFile("in.bam") as source,
        AlignmentFile("out.bam", "wb", template=source) as sink,
        StreamingPileupBuilder(source, tap=sink.write) as builder,
    ):
        pileup = builder.pileup("chr1", 100)
    ```
    """

    def __init__(
        self,
        records: Iterable[AlignedSegment],
        *,
        min_mapq: int = 0,
        min_base_quality: int = DEFAULT_MIN_BASE_QUALITY,
        proper_pairs_only: bool = False,
        include_secondary: bool = False,
        include_supplementary: bool = False,
        include_duplicate: bool = False,
        include_qcfail: bool = False,
        tap: Callable[[AlignedSegment], Any] | None = None,
    ) -> None:
        """Start a builder over coordinate-sorted reads.

        Args:
            records: coordinate-sorted reads, such as an open `AlignmentFile`.
            min_mapq: the lowest mapping quality of a read to pile up.
            min_base_quality: the quality floor of each pileup's filtered views.
            proper_pairs_only: pile up only reads flagged as in a proper pair.
            include_secondary: pile up secondary alignments.
            include_supplementary: pile up supplementary alignments.
            include_duplicate: pile up reads flagged as duplicates.
            include_qcfail: pile up reads flagged as failing quality checks.
            tap: a function given every read once the builder has moved past it.

        Raises:
            ValueError: if the header of the first read does not declare coordinate order.
        """
        self._records: Iterator[AlignedSegment] = iter(records)
        self._next: AlignedSegment | None = next(self._records, None)
        self.header: AlignmentHeader | None = None if self._next is None else self._next.header
        if self.header is not None:
            sort_order = SORT_ORDER.search(str(self.header))
            if sort_order is None or sort_order.group(1) != "coordinate":
                raise ValueError("Records must be coordinate sorted.")
        self.min_mapq: int = min_mapq
        self.min_base_quality: int = min_base_quality
        self.proper_pairs_only: bool = proper_pairs_only
        self.include_secondary: bool = include_secondary
        self.include_supplementary: bool = include_supplementary
        self.include_duplicate: bool = include_duplicate
        self.include_qcfail: bool = include_qcfail
        self.tap: Callable[[AlignedSegment], Any] | None = tap
        self.previous_pileup: Pileup | None = None
        self._waiting: deque[tuple[AlignedSegment, int, int]] = deque()
        self._active: list[Footprint] = []
        self._active_reference_id: int = -1
        self._last_key: tuple[int, int] = (-1, -1)
        self._at: tuple[int, int] | None = None

    def __enter__(self) -> Self:
        """Enter the builder's context."""
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        """Close the builder."""
        self.close()

    def close(self) -> None:
        """Stop, first handing every read not yet handed to `tap` to it, in input order.

        With a `tap`, the rest of the input is read to the end, so an output written by `tap` is
        complete. Without one, no more of the input is read.
        """
        if self.tap is not None:
            while self._waiting:
                self._evict(self._waiting.popleft()[0])
            while self._next is not None:
                self._evict(self._next)
                self._next = next(self._records, None)
        self._waiting.clear()
        self._next = None
        self._active = []

    def accepts(self, record: AlignedSegment) -> bool:
        """Whether a read passes the builder's read filters and has a base on the reference."""
        if record.mapping_quality < self.min_mapq:
            return False
        if self.proper_pairs_only and not record.is_proper_pair:
            return False
        if not self.include_secondary and record.is_secondary:
            return False
        if not self.include_supplementary and record.is_supplementary:
            return False
        if not self.include_duplicate and record.is_duplicate:
            return False
        if not self.include_qcfail and record.is_qcfail:
            return False
        return is_placed(record)

    def pileup(self, contig: str, pos: int) -> Pileup:
        """Advance to a position, at or after the last one, and pile up the reads there.

        Args:
            contig: the name of the contig.
            pos: the 0-based position on the contig.

        Raises:
            ValueError: if the contig is not in the header, the position is negative, or the
                position is before the last one asked for.
        """
        previous = self.previous_pileup
        if (
            previous is not None
            and previous.reference_pos == pos
            and previous.reference_name == contig
        ):
            return previous
        if pos < 0:
            raise ValueError(f"Position must be non-negative, found: {pos}")
        if self.header is None:
            pileup = Pileup(contig, pos, [], self.min_base_quality)
            self.previous_pileup = pileup
            return pileup
        reference_id = self.header.get_tid(contig)
        if reference_id < 0:
            raise ValueError(f"Contig {contig} is not in the header.")
        if self._at is not None and (reference_id, pos) < self._at:
            at = "" if previous is None else f"{previous.reference_name}:{previous.reference_pos}"
            raise ValueError(f"Attempted to advance to {contig}:{pos} from {at}.")
        self._at = (reference_id, pos)
        self._advance(reference_id, pos)
        pileup = Pileup(contig, pos, self._entries(pos), self.min_base_quality)
        self.previous_pileup = pileup
        return pileup

    def columns(self, contig: str, start: int, end: int) -> Iterator[Pileup]:
        """Yield the pileup at every position from `start` to `end`, covered or not.

        Args:
            contig: the name of the contig.
            start: the 0-based first position.
            end: the 0-based position after the last one.

        Raises:
            ValueError: if `end` is before `start`, or as `pileup` does.
        """
        if end < start:
            raise ValueError(f"End {end} is before start {start}.")
        for pos in range(start, end):
            yield self.pileup(contig, pos)

    def _evict(self, record: AlignedSegment) -> None:
        if self.tap is not None:
            self.tap(record)

    def _advance(self, reference_id: int, pos: int) -> None:
        if reference_id != self._active_reference_id:
            self._active = []
            self._active_reference_id = reference_id
        else:
            self._active = [footprint for footprint in self._active if footprint.end > pos]
        target = (reference_id, pos + 1)
        waiting = self._waiting
        tapped = self.tap is not None
        while self._next is not None:
            record = self._next
            record_id = record.reference_id if record.reference_id >= 0 else UNPLACED
            key = (record_id, record.reference_start)
            if key > target:
                break
            if key < self._last_key:
                raise ValueError(f"Records are out of coordinate order at {record.query_name}.")
            self._last_key = key
            if tapped:
                end = record.reference_end
                waiting.append((record, record_id, key[1] + 1 if end is None else end))
            if record_id == reference_id and self.accepts(record):
                footprint = Footprint(record)
                if footprint.end > pos:
                    self._active.append(footprint)
            self._next = next(self._records, None)
        while waiting and (
            waiting[0][1] < reference_id or (waiting[0][1] == reference_id and waiting[0][2] <= pos)
        ):
            self._evict(waiting.popleft()[0])

    def _entries(self, pos: int) -> list[PileupRead]:
        entries: list[PileupRead] = []
        append = entries.append
        for footprint in self._active:
            record = footprint.record
            index = pos - footprint.start
            if index >= 0:
                offset = footprint.offsets[index]
                if offset >= 0:
                    append(PileupRead(record, offset, offset, BASE))
                elif offset == DELETED_AT_END:
                    append(PileupRead(record, None, None, DELETION))
                elif offset != SKIPPED:
                    append(PileupRead(record, None, -offset - 3, DELETION))
            if footprint.insertions:
                insertion = footprint.insertions.get(pos)
                if insertion is not None:
                    append(PileupRead(record, None, None, INSERTION, insertion[0], insertion[1]))
        return entries
