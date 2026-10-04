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

from streampile._footprint import Footprint
from streampile._footprint import is_placed
from streampile._pileup import DEFAULT_MIN_BASE_QUALITY
from streampile._pileup import Pileup
from streampile._pileup import PileupRead
from streampile._pileup import pileup_entries

SORT_ORDER = re.compile(r"^@HD\t.*\bSO:([^\t\n]+)", re.MULTILINE)

UNPLACED: int = sys.maxsize
"""The contig index given to reads with no contig, which sort after every placed read."""


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
        read_filter: Callable[[AlignedSegment], bool] | None = None,
        tap: Callable[[AlignedSegment], Any] | None = None,
    ) -> None:
        """Start a builder over coordinate-sorted reads.

        The header is the `header` of `records` when it has one, such as an `AlignmentFile`, or
        else the first read's, so an empty `AlignmentFile` is checked too. With no header at all,
        e.g. from an empty list, every pileup is empty, and contigs are ordered as they are first
        asked for.

        Args:
            records: coordinate-sorted reads, such as an open `AlignmentFile`.
            min_mapq: the lowest mapping quality of a read to pile up.
            min_base_quality: the quality floor of each pileup's filtered views: 13 by default,
                as in htslib, where `tabulate` counts bases of any quality by default.
            proper_pairs_only: pile up only reads flagged as in a proper pair.
            include_secondary: pile up secondary alignments.
            include_supplementary: pile up supplementary alignments.
            include_duplicate: pile up reads flagged as duplicates.
            include_qcfail: pile up reads flagged as failing quality checks. They are left out
                by default, as htslib leaves them out, although fgbio keeps them.
            read_filter: a function that keeps a read for pileups when it returns True, e.g.
                `lambda read: read.is_proper_pair`, asked only of reads that pass the other
                filters. A read it rejects still goes to `tap`.
            tap: a function given every read once the builder has moved past it.

        Raises:
            ValueError: if the header does not declare coordinate order.
        """
        source_header: object = getattr(records, "header", None)
        self._records: Iterator[AlignedSegment] = iter(records)
        self._next: AlignedSegment | None = next(self._records, None)
        self.header: AlignmentHeader | None = None
        if isinstance(source_header, AlignmentHeader):
            self.header = source_header
        elif self._next is not None:
            self.header = self._next.header
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
        self.read_filter: Callable[[AlignedSegment], bool] | None = read_filter
        self.tap: Callable[[AlignedSegment], Any] | None = tap
        self.previous_pileup: Pileup | None = None
        self._waiting: deque[tuple[AlignedSegment, int, int]] = deque()
        self._active: list[Footprint] = []
        self._active_reference_id: int = -1
        self._last_key: tuple[int, int] = (-1, -1)
        self._at: tuple[int, int] | None = None
        self._contigs_asked: dict[str, int] = {}
        self._closed: bool = False

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
        tap = self.tap
        if tap is not None:
            while self._waiting:
                tap(self._waiting.popleft()[0])
            while self._next is not None:
                tap(self._next)
                self._next = next(self._records, None)
        self._waiting.clear()
        self._next = None
        self._active = []
        self._closed = True

    def accepts(self, record: AlignedSegment) -> bool:
        """Whether a read passes the built-in filters, is placed, and then passes `read_filter`."""
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
        if not is_placed(record):
            return False
        return self.read_filter is None or self.read_filter(record)

    def pileup(self, contig: str, pos: int) -> Pileup:
        """Advance to a position, at or after the last one, and pile up the reads there.

        Args:
            contig: the name of the contig.
            pos: the 0-based position on the contig.

        Raises:
            ValueError: if the builder is closed, the contig is not in the header, the position is
                negative, or the position is before the last one asked for.
        """
        if self._closed:
            raise ValueError("The builder is closed.")
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
            reference_id = self._contigs_asked.setdefault(contig, len(self._contigs_asked))
        else:
            reference_id = self.header.get_tid(contig)
            if reference_id < 0:
                raise ValueError(f"Contig {contig} is not in the header.")
        if self._at is not None and (reference_id, pos) < self._at:
            at = "" if previous is None else f"{previous.reference_name}:{previous.reference_pos}"
            raise ValueError(f"Attempted to advance to {contig}:{pos} from {at}.")
        self._at = (reference_id, pos)
        if self.header is None:
            entries: list[PileupRead] = []
        else:
            self._advance(reference_id, pos)
            entries = pileup_entries(self._active, pos)
        pileup = Pileup(contig, pos, entries, self.min_base_quality)
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

    def _advance(self, reference_id: int, pos: int) -> None:
        if reference_id != self._active_reference_id:
            self._active = []
            self._active_reference_id = reference_id
        else:
            self._active = [footprint for footprint in self._active if footprint.end > pos]
        target = (reference_id, pos + 1)
        waiting = self._waiting
        tap = self.tap
        while self._next is not None:
            record = self._next
            record_id = record.reference_id if record.reference_id >= 0 else UNPLACED
            key = (record_id, record.reference_start)
            if key > target:
                break
            if key < self._last_key:
                raise ValueError(f"Records are out of coordinate order at {record.query_name}.")
            self._last_key = key
            if tap is not None:
                end = record.reference_end
                waiting.append((record, record_id, key[1] + 1 if end is None else end))
            if record_id == reference_id and self.accepts(record):
                footprint = Footprint(record)
                if footprint.end > pos:
                    self._active.append(footprint)
            self._next = next(self._records, None)
        if tap is not None:
            while waiting and (
                waiting[0][1] < reference_id
                or (waiting[0][1] == reference_id and waiting[0][2] <= pos)
            ):
                tap(waiting.popleft()[0])
