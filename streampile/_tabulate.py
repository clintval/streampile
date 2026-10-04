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

from streampile._footprint import query_qualities
from streampile._pileup import DEFAULT_EXCLUDE_FLAGS
from streampile._table import TabulatedBase

ACGT: Final[frozenset[str]] = frozenset("ACGT")

REFERENCE_PADDING: Final[int] = 10_000

CHUNK: Final[int] = 1_000_000
"""The most bases of a territory span counted at once, which bounds memory on long spans."""


@final
class _Reference:
    """Two sliding windows over one contig of an indexed FASTA, upper-cased.

    A read with a long reference skip reads two far-apart stretches, so the window used last is
    kept beside the other, and neither is fetched again for the next read.
    """

    def __init__(self, fasta: FastaFile, contig: str) -> None:
        self.fasta: FastaFile = fasta
        self.contig: str = contig
        self.length: int = fasta.get_reference_length(contig)
        self.windows: list[tuple[int, int, str]] = []

    def get(self, start: int, end: int) -> str:
        start, end = max(start, 0), min(end, self.length)
        for index, (low, high, sequence) in enumerate(self.windows):
            if low <= start and end <= high:
                if index:
                    self.windows.reverse()
                return sequence[start - low : end - low]
        low = max(start - REFERENCE_PADDING, 0)
        high = min(max(end, start) + 10 * REFERENCE_PADDING, self.length)
        sequence = self.fasta.fetch(self.contig, low, high).upper()
        self.windows = [(low, high, sequence), *self.windows[:1]]
        return sequence[start - low : end - low]


@final
class _Ledger:
    """The counts of one territory span while reads are added to it, by strand."""

    __slots__ = ("alleles", "depth", "end", "no_calls", "ref_fwd", "ref_rev", "start")

    def __init__(self, start: int, end: int) -> None:
        self.start: int = start
        self.end: int = end
        self.depth: list[int] = [0] * (end - start + 1)
        self.ref_fwd: list[int] = [0] * (end - start + 1)
        self.ref_rev: list[int] = [0] * (end - start + 1)
        self.no_calls: list[int] = [0] * (end - start)
        self.alleles: dict[int, dict[tuple[str, str], list[int]]] = {}

    def add(self, start: int, end: int, *, ref: bool, reverse: bool) -> None:
        start, end = max(start, self.start), min(end, self.end)
        if start < end:
            self.depth[start - self.start] += 1
            self.depth[end - self.start] -= 1
            if ref:
                counts = self.ref_rev if reverse else self.ref_fwd
                counts[start - self.start] += 1
                counts[end - self.start] -= 1

    def add_no_call(self, pos: int) -> None:
        if self.start <= pos < self.end:
            self.no_calls[pos - self.start] += 1

    def add_allele(self, pos: int, key: tuple[str, str], *, reverse: bool) -> None:
        if self.start <= pos < self.end:
            self.alleles.setdefault(pos, {}).setdefault(key, [0, 0])[reverse] += 1

    def bases(self, contig: str, reference: _Reference) -> Iterator[TabulatedBase]:
        bases = reference.get(self.start, self.end)
        depth = ref_fwd = ref_rev = 0
        for offset, base in enumerate(bases):
            depth += self.depth[offset]
            ref_fwd += self.ref_fwd[offset]
            ref_rev += self.ref_rev[offset]
            counts = self.alleles.get(self.start + offset)
            if counts is None:
                yield TabulatedBase(
                    contig=contig,
                    pos=self.start + offset + 1,
                    ref=base,
                    depth=depth,
                    no_calls=self.no_calls[offset],
                    ref_reads=ref_fwd + ref_rev,
                    ref_fwd=ref_fwd,
                    ref_rev=ref_rev,
                )
                continue
            alleles = sorted(counts.items(), key=_by_reads)
            yield TabulatedBase(
                contig=contig,
                pos=self.start + offset + 1,
                ref=base,
                depth=depth,
                no_calls=self.no_calls[offset],
                ref_reads=ref_fwd + ref_rev,
                ref_fwd=ref_fwd,
                ref_rev=ref_rev,
                alt_refs=tuple(allele[0] for allele, _ in alleles),
                alts=tuple(allele[1] for allele, _ in alleles),
                alt_reads=tuple(fwd + rev for _, (fwd, rev) in alleles),
                alt_fwd=tuple(fwd for _, (fwd, _) in alleles),
                alt_rev=tuple(rev for _, (_, rev) in alleles),
            )


def _by_reads(item: tuple[tuple[str, str], list[int]]) -> tuple[int, str, str]:
    return -sum(item[1]), item[0][0], item[0][1]


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
    """Adjacent differences of one read, and whether a matching aligned base borders each end.

    `floor` is the lowest position the run's allele may be left-aligned to: the end of the read's
    previous difference or reference skip, or else the read's first aligned base.
    """

    events: list[_Event]
    anchored: bool
    closed: bool
    floor: int


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
    """Trim and left-align an allele as `bcftools norm` does.

    Returns `None` for an allele that changes nothing, or that would move past `floor`.

    Args:
        pos: the 0-based position of the allele's first reference base.
        ref: the reference bases.
        alt: the alternate bases.
        reference: the bases of the contig from a 0-based start to an end.
        floor: the lowest position the allele may move to.
    """
    if ref == alt:
        return None
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

    def __init__(self, start: int) -> None:
        self.runs: list[_Run] = []
        self.events: list[_Event] = []
        self.anchored: bool = False
        self.floor: int = start

    def add(self, event: _Event) -> None:
        self.events.append(event)

    def border(self, *, aligned: bool) -> None:
        """End any run at a matching aligned base, or at a clip, skip, or end of the read."""
        if self.events:
            self.runs.append(_Run(self.events, self.anchored, aligned, self.floor))
            self.floor = self.events[-1].ref_end
            self.events = []
        self.anchored = aligned


def _runs(
    segments: list[Segment], sequence: str, reference: "_Reference", start: int
) -> list[_Run]:
    """The runs of adjacent differences from the reference of one read, in alignment order."""
    runs = _Runs(start)
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
            if operator == CREF_SKIP:
                runs.floor = ref_pos + length
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
    A run that spells the reference, such as an insertion and a deletion of the same base, is no
    allele, and its read is a reference read across it.

    A read is counted for an allele only when none of its bases in the allele is an `N`, and its
    mismatched and inserted bases are at the base-quality floor, as is the base before an insertion
    that opens the allele and the base after a deletion that closes it: as in htslib, an insertion
    is judged with the base it follows, and a deletion by the read's next base, not its anchor.
    Otherwise the read is not informative at any base of that allele. A read is never counted for an
    allele that starts with an indel with no aligned base before it, ends with a deletion with no
    aligned base after it, or holds a reference base other than A, C, G, or T, before or after
    left-alignment. Nor is it counted for an indel that would be left-aligned past the read's
    previous difference or reference skip, or past its first aligned base, so one read is never
    counted for two alleles at one base; it is then not informative from there to the end of the
    indel. A read with no stored qualities (QUAL `*`) has quality 255 at every base, as in htslib,
    so it passes every floor, and a read with no stored bases (SEQ `*`) is not counted at all.

    A read counted for an allele is informative at every base the allele spans, a deletion's
    deleted bases included, and at every base an indel is left-aligned across. Elsewhere, a read
    is informative at a base where it holds an aligned base at the floor that is not an `N`, over
    a reference base of A, C, G, or T. Reference skips, the CIGAR `N` operator, and clips cover
    nothing. A read's `N` bases, its no-calls, are counted apart, whatever their quality.
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
        self._low_quality: bytes = bytes(int(quality < min_base_quality) for quality in range(256))
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
        """Iterate over every base of a territory, covered or not.

        Bases come in the order of the alignment header's contigs, then by position. Every contig
        of the territory is checked against the alignment header and the reference when this is
        called, before any base is read, so a refused territory leaves nothing half written.

        Args:
            alignments: an indexed, coordinate-sorted BAM or CRAM.
            territory: the bases to tabulate.

        Raises:
            ValueError: if a contig of the territory is not in the alignment header or the
                reference, or has another length in each, or a span runs past the end of its
                contig.
        """
        contigs = _by_header(alignments, self.reference, territory)
        return (
            base
            for contig, spans in contigs
            for base in self._tabulate_contig(alignments, contig, spans)
        )

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
        reference = self._reference(record.reference_name)
        counted, dropped, _ = self._alleles(record, _segments(record), reference)
        return counted, dropped

    def _alleles(
        self, record: AlignedSegment, segments: list[Segment], reference: _Reference
    ) -> tuple[list[Allele], list[tuple[int, int]], list[tuple[int, int]]]:
        """The alleles of a read, the stretches it is not counted for, and those it matches."""
        sequence: str = record.query_sequence or ""
        qualities = query_qualities(record) if self.min_base_quality > 0 else None
        counted: list[Allele] = []
        dropped: list[tuple[int, int]] = []
        matched: list[tuple[int, int]] = []
        for run in _runs(segments, sequence, reference, record.reference_start):
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
            opens_with_insertion = first.is_indel and first.ref_end == first.ref_start
            closes_with_deletion = last.is_indel and last.ref_end > last.ref_start
            judged_from = query_start if opens_with_insertion else first.query_start
            judged_to = query_end + 1 if closes_with_deletion else query_end
            if (
                (closes_with_deletion and not run.closed)
                or (
                    qualities is not None
                    and min(qualities[judged_from:judged_to]) < self.min_base_quality
                )
                or "N" in alt
                or not set(ref) <= ACGT
            ):
                dropped.append((ref_start, ref_end))
                continue
            if ref == alt:
                matched.append((ref_start, ref_end))
                continue
            normalized = normalize(ref_start, ref, alt, reference.get, run.floor)
            if normalized is None:
                dropped.append((run.floor, ref_end))
                continue
            pos, ref, alt = normalized
            start, end = min(ref_start, pos), max(ref_end, pos + len(ref))
            if reference.get(start, ref_start).strip("ACGT"):
                dropped.append((start, end))
                continue
            counted.append(Allele(pos=pos, ref=ref, alt=alt, start=start, end=end))
        return counted, dropped, matched

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
        reverse = record.is_reverse
        counted, dropped, matched = self._alleles(record, segments, reference)
        uninformative, no_calls = self._uninformative(record, segments, reference)
        excluded: set[int] = {*uninformative, *no_calls}
        for pos in no_calls:
            for ledger in touched:
                ledger.add_no_call(pos)
        for start, end in dropped:
            excluded.update(range(start, end))
        for allele in counted:
            after = allele.pos + len(allele.ref)
            _add(touched, allele.start, allele.pos, ref=True, reverse=reverse)
            _add(touched, allele.pos, after, ref=False, reverse=reverse)
            _add(touched, after, allele.end, ref=True, reverse=reverse)
            excluded.update(range(allele.start, allele.end))
            for ledger in touched:
                ledger.add_allele(allele.pos, (allele.ref, allele.alt), reverse=reverse)
        for start, end in matched:
            _add(touched, start, end, ref=True, reverse=reverse)
            excluded.update(range(start, end))
        for operator, block_start, _, length in segments:
            if operator == CMATCH:
                cursor = block_start
                for pos in sorted(p for p in excluded if block_start <= p < block_start + length):
                    _add(touched, cursor, pos, ref=True, reverse=reverse)
                    cursor = pos + 1
                _add(touched, cursor, block_start + length, ref=True, reverse=reverse)

    def _uninformative(
        self, record: AlignedSegment, segments: list[Segment], reference: _Reference
    ) -> tuple[list[int], list[int]]:
        """The aligned positions of a read under the floor or over a non-ACGT base, and its Ns."""
        sequence: str = record.query_sequence or ""
        qualities = query_qualities(record)
        low = None
        if qualities is not None and self.min_base_quality > 0:
            low = bytes(qualities).translate(self._low_quality)
        uninformative: list[int] = []
        no_calls: list[int] = []
        for operator, ref_pos, query, length in segments:
            if operator != CMATCH:
                continue
            end = query + length
            at = sequence.find("N", query, end)
            while at >= 0:
                no_calls.append(ref_pos + at - query)
                at = sequence.find("N", at + 1, end)
            if low is not None:
                at = low.find(1, query, end)
                while at >= 0:
                    uninformative.append(ref_pos + at - query)
                    at = low.find(1, at + 1, end)
            bases = reference.get(ref_pos, ref_pos + length)
            if bases.strip("ACGT"):
                uninformative.extend(
                    ref_pos + offset for offset, base in enumerate(bases) if base not in ACGT
                )
        return uninformative, no_calls


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


def _add(ledgers: list[_Ledger], start: int, end: int, *, ref: bool, reverse: bool) -> None:
    for ledger in ledgers:
        ledger.add(start, end, ref=ref, reverse=reverse)


def _contig_length(alignments: AlignmentFile, reference: FastaFile, contig: str) -> int:
    """The length of a contig, checked against the alignment header and the reference.

    Raises:
        ValueError: if the contig is not in the alignment header or the reference, or has another
            length in each.
    """
    if alignments.get_tid(contig) < 0:
        raise ValueError(f"Contig {contig} is not in the alignment header.")
    if contig not in reference.references:
        raise ValueError(f"Contig {contig} is not in the reference.")
    length = alignments.get_reference_length(contig)
    if reference.get_reference_length(contig) != length:
        raise ValueError(
            f"Contig {contig} has {reference.get_reference_length(contig)} bases in the reference"
            + f" but {length} in the alignment header."
        )
    return length


def _by_header(
    alignments: AlignmentFile, reference: FastaFile, territory: Territory
) -> list[tuple[str, list[tuple[int, int]]]]:
    """The spans of a territory, grouped by contig in the order of the alignment header.

    Raises:
        ValueError: if a contig of the territory is not in the alignment header or the
            reference, or has another length in each, or a span runs past the end of its contig.
    """
    lengths: dict[str, int] = {}
    by_contig: dict[str, list[tuple[int, int]]] = {}
    for span in territory:
        contig = span.refname
        length = lengths.get(contig)
        if length is None:
            length = lengths[contig] = _contig_length(alignments, reference, contig)
        if span.end > length:
            raise ValueError(
                f"Span {contig}:{span.start}-{span.end} runs past the end of {contig},"
                + f" which has {length} bases."
            )
        by_contig.setdefault(contig, []).extend(
            (start, min(start + CHUNK, span.end)) for start in range(span.start, span.end, CHUNK)
        )
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
    """Iterate over every base of a territory with its depth and the reads of each allele.

    See `Tabulator` for which reads count where. Bases come in the order of the alignment
    header's contigs, then by position. The territory is checked when this is called, as in
    `Tabulator.tabulate`.

    Args:
        alignments: an indexed, coordinate-sorted BAM or CRAM.
        reference: the indexed reference the reads are aligned to.
        territory: the bases to tabulate, e.g. `Territory(BedReader.from_path[Bed3N](path))`.
        min_base_quality: the lowest base quality of an informative base: 0 by default, where a
            `StreamingPileupBuilder` filters at 13 by default, as htslib does.
        min_mapping_quality: the lowest mapping quality of a counted read.
        exclude_flags: reads with any of these SAM flags are not counted: by default, secondary,
            QC-fail, duplicate, and supplementary reads.

    Raises:
        ValueError: if a contig of the territory is not in the alignment header or the
            reference, or has another length in each, or a span runs past the end of its contig.
    """
    tabulator = Tabulator(
        reference,
        min_base_quality=min_base_quality,
        min_mapping_quality=min_mapping_quality,
        exclude_flags=exclude_flags,
    )
    return tabulator.tabulate(alignments, territory)
