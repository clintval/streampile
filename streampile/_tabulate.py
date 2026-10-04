from bisect import bisect_right
from collections.abc import Callable
from collections.abc import Iterator
from dataclasses import dataclass
from typing import Final
from typing import final

from bedspec import Territory
from pysam import CDEL
from pysam import CDIFF
from pysam import CEQUAL
from pysam import CINS
from pysam import CMATCH
from pysam import CREF_SKIP
from pysam import CSOFT_CLIP
from pysam import AlignedSegment
from pysam import AlignmentFile
from pysam import FastaFile

from streampile._table import TabulatedBase

DEFAULT_EXCLUDE_FLAGS: Final[int] = 0xF00
"""Secondary, QC-failed, duplicate, and supplementary reads, which are left out by default."""

ACGT: Final[frozenset[str]] = frozenset("ACGT")

REFERENCE_PADDING: Final[int] = 10_000


@final
class _Reference:
    """A sliding window over one contig of an indexed FASTA, upper-cased."""

    def __init__(self, fasta: FastaFile, contig: str) -> None:
        self.fasta: FastaFile = fasta
        self.contig: str = contig
        self.length: int = fasta.get_reference_length(contig)
        self.start: int = 0
        self.end: int = 0
        self.sequence: str = ""

    def get(self, start: int, end: int) -> str:
        start, end = max(start, 0), min(end, self.length)
        if start < self.start or end > self.end:
            self.start = max(start - REFERENCE_PADDING, 0)
            self.end = min(max(end, start) + 10 * REFERENCE_PADDING, self.length)
            self.sequence = self.fasta.fetch(self.contig, self.start, self.end).upper()
        return self.sequence[start - self.start : end - self.start]


@final
class _Ledger:
    """The counts of one territory span while reads are added to it."""

    __slots__ = ("alleles", "depth", "end", "ref_reads", "start")

    def __init__(self, start: int, end: int) -> None:
        self.start: int = start
        self.end: int = end
        self.depth: list[int] = [0] * (end - start + 1)
        self.ref_reads: list[int] = [0] * (end - start + 1)
        self.alleles: dict[int, dict[tuple[str, str], int]] = {}

    def add(self, start: int, end: int, *, ref: bool) -> None:
        start, end = max(start, self.start), min(end, self.end)
        if start < end:
            self.depth[start - self.start] += 1
            self.depth[end - self.start] -= 1
            if ref:
                self.ref_reads[start - self.start] += 1
                self.ref_reads[end - self.start] -= 1

    def add_allele(self, pos: int, key: tuple[str, str]) -> None:
        if self.start <= pos < self.end:
            counts = self.alleles.setdefault(pos, {})
            counts[key] = counts.get(key, 0) + 1

    def bases(self, contig: str, reference: _Reference) -> Iterator[TabulatedBase]:
        bases = reference.get(self.start, self.end)
        depth = 0
        ref_reads = 0
        for offset, base in enumerate(bases):
            depth += self.depth[offset]
            ref_reads += self.ref_reads[offset]
            counts = self.alleles.get(self.start + offset)
            if counts is None:
                yield TabulatedBase(
                    contig=contig,
                    pos=self.start + offset + 1,
                    ref=base,
                    depth=depth,
                    ref_reads=ref_reads,
                )
                continue
            alleles = sorted(counts.items(), key=_by_reads)
            yield TabulatedBase(
                contig=contig,
                pos=self.start + offset + 1,
                ref=base,
                depth=depth,
                ref_reads=ref_reads,
                alt_refs=tuple(allele[0] for allele, _ in alleles),
                alts=tuple(allele[1] for allele, _ in alleles),
                alt_reads=tuple(reads for _, reads in alleles),
            )


def _by_reads(item: tuple[tuple[str, str], int]) -> tuple[int, str, str]:
    return -item[1], item[0][0], item[0][1]


@dataclass(frozen=True, slots=True)
class _Event:
    """A difference from the reference within one read, in reference and query coordinates."""

    ref_start: int
    ref_end: int
    query_start: int
    query_end: int
    is_indel: bool


@dataclass(frozen=True, slots=True)
class _Run:
    """Adjacent differences of one read, and whether a matching aligned base borders each end."""

    events: list[_Event]
    anchored: bool
    closed: bool


@dataclass(frozen=True, slots=True)
class Allele:
    """A normalized VCF allele seen in one read, and the reference stretch it accounts for.

    Attributes:
        pos: the 0-based position of the allele's first reference base.
        ref: the reference bases.
        alt: the read's bases.
        start: the 0-based first position the read is informative at because of the allele.
        end: the 0-based position after the last one.
    """

    pos: int
    ref: str
    alt: str
    start: int
    end: int

    @property
    def key(self) -> str:
        """The allele as `REF>ALT`."""
        return f"{self.ref}>{self.alt}"


def normalize(
    pos: int, ref: str, alt: str, reference: Callable[[int, int], str], floor: int = 0
) -> tuple[int, str, str] | None:
    """Trim and left-align an allele as `bcftools norm` does, or `None` if it would pass `floor`.

    Args:
        pos: the 0-based position of the allele's first reference base.
        ref: the reference bases.
        alt: the alternate bases.
        reference: the bases of the contig from a 0-based start to an end.
        floor: the lowest position the allele may move to.
    """
    while ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
        if not ref or not alt:
            if pos - 1 < floor:
                return None
            pos -= 1
            base = reference(pos, pos + 1)
            ref, alt = base + ref, base + alt
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
        pos += 1
    return pos, ref, alt


Segment = tuple[int, int, int, int]
"""One CIGAR operator of a read: the operator, its reference start, its query start, and length."""


def _segments(record: AlignedSegment) -> list[Segment]:
    segments: list[Segment] = []
    ref_pos: int = record.reference_start
    query = 0
    for operator, length in record.cigartuples or ():
        if operator == CMATCH or operator == CEQUAL or operator == CDIFF:
            segments.append((CMATCH, ref_pos, query, length))
            ref_pos += length
            query += length
        elif operator == CINS:
            segments.append((CINS, ref_pos, query, length))
            query += length
        elif operator == CDEL or operator == CREF_SKIP:
            segments.append((operator, ref_pos, query, length))
            ref_pos += length
        elif operator == CSOFT_CLIP:
            segments.append((CSOFT_CLIP, ref_pos, query, length))
            query += length
    return segments


@final
class _Runs:
    """Collects the runs of one read as its aligned bases and indels are walked in order."""

    def __init__(self) -> None:
        self.runs: list[_Run] = []
        self.events: list[_Event] = []
        self.anchored: bool = False

    def add(self, event: _Event) -> None:
        self.events.append(event)

    def border(self, *, aligned: bool) -> None:
        """End any run at a matching aligned base, or at a clip, skip, or end of the read."""
        if self.events:
            self.runs.append(_Run(self.events, self.anchored, aligned))
            self.events = []
        self.anchored = aligned


def _runs(segments: list[Segment], sequence: str, reference: "_Reference") -> list[_Run]:
    """The runs of adjacent differences from the reference of one read, in alignment order."""
    runs = _Runs()
    for operator, ref_pos, query, length in segments:
        if operator == CMATCH:
            bases = reference.get(ref_pos, ref_pos + length)
            read = sequence[query : query + length]
            previous = -1
            if read != bases:
                for offset, (base, ref_base) in enumerate(zip(read, bases, strict=False)):
                    if base != ref_base:
                        if offset > previous + 1:
                            runs.border(aligned=True)
                        at, query_at = ref_pos + offset, query + offset
                        runs.add(_Event(at, at + 1, query_at, query_at + 1, False))
                        previous = offset
            if previous < len(bases) - 1:
                runs.border(aligned=True)
            if len(bases) < length:
                runs.border(aligned=False)
        elif operator == CINS:
            runs.add(_Event(ref_pos, ref_pos, query, query + length, True))
        elif operator == CDEL:
            runs.add(_Event(ref_pos, ref_pos + length, query, query, True))
        else:
            runs.border(aligned=False)
    runs.border(aligned=False)
    return runs.runs


class Tabulator:
    """Count the reads of every allele at every base of a territory.

    Reads are filtered by flag and mapping quality, then each read is walked once. Within a read,
    every run of adjacent differences from the reference, with no matching aligned base between
    them, is one allele: one mismatch is an SNV, adjacent mismatches are an MNV, and an indel,
    alone or with mismatches next to it, is an indel or a complex allele. An allele that starts
    with an indel is anchored on the aligned base before it, as in VCF. Alleles are then trimmed
    and left-aligned, like `bcftools norm`, so they can be matched to other VCF alleles exactly.

    A read is counted for an allele only when every read base of the allele, the anchor base
    included, is at the base-quality floor and is not an `N`. Otherwise the read is not
    informative at any base of that allele. A read is never counted for an allele that starts
    with an indel with no aligned base before it, ends with a deletion with no aligned base after
    it, would be left-aligned to before the read's first aligned base, or holds a reference base
    other than A, C, G, or T. A read with no stored qualities (QUAL `*`) has quality 255 at every
    base, as in htslib, so it passes every floor, and a read with no stored bases (SEQ `*`) is
    not counted at all.

    A read counted for an allele is informative at every base the allele spans, a deletion's
    deleted bases included, and at every base an indel is left-aligned across. Elsewhere, a read
    is informative at a base where it holds an aligned base at the floor that is not an `N`, over
    a reference base of A, C, G, or T. Reference skips (`N` operators) and clips cover nothing.
    """

    def __init__(
        self,
        reference: FastaFile,
        *,
        min_base_quality: int = 0,
        min_mapping_quality: int = 0,
        exclude_flags: int = DEFAULT_EXCLUDE_FLAGS,
    ) -> None:
        """Set the reference and the read and base filters.

        Args:
            reference: the indexed reference the reads are aligned to.
            min_base_quality: the lowest base quality of an informative base: 0 by default, where
                a `StreamingPileupBuilder` filters at 13 by default, as htslib does.
            min_mapping_quality: the lowest mapping quality of a counted read.
            exclude_flags: reads with any of these SAM flags are not counted: by default,
                secondary, QC-fail, duplicate, and supplementary reads. Unmapped reads never are.
        """
        self.reference: FastaFile = reference
        self.min_base_quality: int = min_base_quality
        self.min_mapping_quality: int = min_mapping_quality
        self.exclude_flags: int = exclude_flags
        self._floor: str = chr(min_base_quality + 33)
        self._low_quality: bytes = bytes(int(code - 33 < min_base_quality) for code in range(256))
        self._contigs: dict[str, _Reference] = {}

    def accepts(self, record: AlignedSegment) -> bool:
        """Whether a read passes the flag and mapping-quality filters and has bases to count."""
        return (
            not record.is_unmapped
            and not record.flag & self.exclude_flags
            and record.mapping_quality >= self.min_mapping_quality
            and record.query_sequence is not None
            and bool(record.cigartuples)
        )

    def tabulate(self, alignments: AlignmentFile, territory: Territory) -> Iterator[TabulatedBase]:
        """Yield every base of a territory, covered or not.

        Bases come in the order of the alignment header's contigs, then by position.

        Args:
            alignments: an indexed, coordinate-sorted alignment file.
            territory: the bases to tabulate.

        Raises:
            ValueError: if a contig of the territory is not in the alignment header.
        """
        for contig, spans in _by_header(alignments, territory):
            yield from self._tabulate_contig(alignments, contig, spans)

    def alleles(self, record: AlignedSegment) -> tuple[list[Allele], list[tuple[int, int]]]:
        """The alleles a read is counted for, and the stretches of those it is not counted for.

        Args:
            record: a mapped read with a sequence.

        Returns:
            The read's counted alleles, and the 0-based half-open reference stretches of the
            alleles it is not counted for, where it is not informative.

        Raises:
            ValueError: if the read is not mapped.
        """
        if record.reference_name is None:
            raise ValueError(f"Read {record.query_name} is not mapped.")
        return self._alleles(record, _segments(record), self._reference(record.reference_name))

    def _alleles(
        self, record: AlignedSegment, segments: list[Segment], reference: _Reference
    ) -> tuple[list[Allele], list[tuple[int, int]]]:
        sequence: str = record.query_sequence or ""
        qualities: str | None = record.query_qualities_str
        counted: list[Allele] = []
        dropped: list[tuple[int, int]] = []
        for run in _runs(segments, sequence, reference):
            first, last = run.events[0], run.events[-1]
            ref_start, query_start = first.ref_start, first.query_start
            ref_end, query_end = last.ref_end, last.query_end
            if first.is_indel:
                if not run.anchored:
                    dropped.append((ref_start, ref_end))
                    continue
                ref_start, query_start = ref_start - 1, query_start - 1
            ref = reference.get(ref_start, ref_end)
            alt = sequence[query_start:query_end]
            if (
                (last.is_indel and last.ref_end > last.ref_start and not run.closed)
                or (qualities is not None and min(qualities[query_start:query_end]) < self._floor)
                or "N" in alt
                or not set(ref) <= ACGT
                or (
                    normalized := normalize(
                        ref_start, ref, alt, reference.get, record.reference_start
                    )
                )
                is None
            ):
                dropped.append((ref_start, ref_end))
                continue
            pos, ref, alt = normalized
            counted.append(
                Allele(
                    pos=pos,
                    ref=ref,
                    alt=alt,
                    start=min(ref_start, pos),
                    end=max(ref_end, pos + len(ref)),
                )
            )
        return counted, dropped

    def _reference(self, contig: str) -> _Reference:
        reference = self._contigs.get(contig)
        if reference is None:
            reference = self._contigs[contig] = _Reference(self.reference, contig)
        return reference

    def _tabulate_contig(
        self, alignments: AlignmentFile, contig: str, spans: list[tuple[int, int]]
    ) -> Iterator[TabulatedBase]:
        reference = self._reference(contig)
        starts = [start for start, _ in spans]
        ledgers: dict[int, _Ledger] = {}
        fetched_to = -1
        for index, (start, end) in enumerate(spans):
            for record in alignments.fetch(contig, start, end):
                if record.reference_start >= fetched_to and self.accepts(record):
                    self._count(record, reference, spans, starts, ledgers)
            fetched_to = end
            ledger = ledgers.pop(index, None) or _Ledger(start, end)
            yield from ledger.bases(contig, reference)

    def _count(
        self,
        record: AlignedSegment,
        reference: _Reference,
        spans: list[tuple[int, int]],
        starts: list[int],
        ledgers: dict[int, _Ledger],
    ) -> None:
        touched = _touched(record, spans, starts, ledgers)
        if not touched:
            return
        segments = _segments(record)
        counted, dropped = self._alleles(record, segments, reference)
        excluded: set[int] = set(self._uninformative(record, segments, reference))
        for start, end in dropped:
            excluded.update(range(start, end))
        claimed: set[int] = set()
        for allele in counted:
            for pos in range(allele.start, allele.end):
                if pos not in claimed:
                    claimed.add(pos)
                    is_ref = not allele.pos <= pos < allele.pos + len(allele.ref)
                    _add(touched, pos, pos + 1, ref=is_ref)
            for ledger in touched:
                ledger.add_allele(allele.pos, (allele.ref, allele.alt))
        excluded |= claimed
        for operator, block_start, _, length in segments:
            if operator == CMATCH:
                cursor = block_start
                for pos in sorted(p for p in excluded if block_start <= p < block_start + length):
                    _add(touched, cursor, pos, ref=True)
                    cursor = pos + 1
                _add(touched, cursor, block_start + length, ref=True)

    def _uninformative(
        self, record: AlignedSegment, segments: list[Segment], reference: _Reference
    ) -> Iterator[int]:
        """The aligned positions of a read below the floor, at a read N, or over a non-ACGT base."""
        sequence: str = record.query_sequence or ""
        qualities: str | None = record.query_qualities_str
        low: list[int] = []
        if qualities is not None and self.min_base_quality > 0:
            mask = qualities.encode().translate(self._low_quality)
            at = mask.find(1)
            while at >= 0:
                low.append(at)
                at = mask.find(1, at + 1)
        at = sequence.find("N")
        while at >= 0:
            low.append(at)
            at = sequence.find("N", at + 1)
        blocks = [segment for segment in segments if segment[0] == CMATCH]
        for offset in low:
            for _, ref_pos, query, length in blocks:
                if query <= offset < query + length:
                    yield ref_pos + offset - query
                    break
        for _, ref_pos, _, length in blocks:
            bases = reference.get(ref_pos, ref_pos + length)
            if not set(bases) <= ACGT:
                for offset, base in enumerate(bases):
                    if base not in ACGT:
                        yield ref_pos + offset


def _touched(
    record: AlignedSegment,
    spans: list[tuple[int, int]],
    starts: list[int],
    ledgers: dict[int, _Ledger],
) -> list[_Ledger]:
    """The ledgers of the spans a read overlaps, opened as needed."""
    read_start: int = record.reference_start
    read_end: int = record.reference_end or read_start
    touched: list[_Ledger] = []
    index = max(bisect_right(starts, read_start) - 1, 0)
    while index < len(spans) and spans[index][0] < read_end:
        if spans[index][1] > read_start:
            ledger = ledgers.get(index)
            if ledger is None:
                ledger = ledgers[index] = _Ledger(*spans[index])
            touched.append(ledger)
        index += 1
    return touched


def _add(ledgers: list[_Ledger], start: int, end: int, *, ref: bool) -> None:
    for ledger in ledgers:
        ledger.add(start, end, ref=ref)


def _by_header(
    alignments: AlignmentFile, territory: Territory
) -> list[tuple[str, list[tuple[int, int]]]]:
    """The spans of a territory, grouped by contig in the order of the alignment header.

    Raises:
        ValueError: if a contig of the territory is not in the alignment header.
    """
    by_contig: dict[str, list[tuple[int, int]]] = {}
    for span in territory:
        if span.refname not in by_contig and alignments.get_tid(span.refname) < 0:
            raise ValueError(f"Contig {span.refname} is not in the alignment header.")
        by_contig.setdefault(span.refname, []).append((span.start, span.end))
    return sorted(by_contig.items(), key=lambda contig: alignments.get_tid(contig[0]))


def tabulate(
    alignments: AlignmentFile,
    reference: FastaFile,
    territory: Territory,
    *,
    min_base_quality: int = 0,
    min_mapping_quality: int = 0,
    exclude_flags: int = DEFAULT_EXCLUDE_FLAGS,
) -> Iterator[TabulatedBase]:
    """Yield every base of a territory with its depth and the reads of each allele.

    See `Tabulator` for which reads count where. Bases come in the order of the alignment
    header's contigs, then by position.

    Args:
        alignments: an indexed, coordinate-sorted alignment file.
        reference: the indexed reference the reads are aligned to.
        territory: the bases to tabulate, e.g. `Territory(BedReader.from_path[Bed3N](path))`.
        min_base_quality: the lowest base quality of an informative base: 0 by default, where a
            `StreamingPileupBuilder` filters at 13 by default, as htslib does.
        min_mapping_quality: the lowest mapping quality of a counted read.
        exclude_flags: reads with any of these SAM flags are not counted: by default, secondary,
            QC-fail, duplicate, and supplementary reads.
    """
    tabulator = Tabulator(
        reference,
        min_base_quality=min_base_quality,
        min_mapping_quality=min_mapping_quality,
        exclude_flags=exclude_flags,
    )
    yield from tabulator.tabulate(alignments, territory)
