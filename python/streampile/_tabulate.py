from collections.abc import Callable
from collections.abc import Iterator
from typing import Final

from bedspec import Territory
from pysam import AlignedSegment
from pysam import AlignmentFile
from pysam import FastaFile

from streampile import _native
from streampile._frozen import frozen
from streampile._pileup import DEFAULT_EXCLUDE_FLAGS
from streampile._table import TabulatedBase

CHUNK: Final[int] = 1_000_000
"""The most bases of a territory span counted at once, which bounds memory on long spans."""


@frozen(slots=True)
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
    return _native.normalize_allele(pos, ref, alt, reference, floor)


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
    that opens the allele and the base after a deletion that closes it: as in pysam's `pileup()` and
    `samtools mpileup`, an insertion is judged with the base it follows, and a deletion by the
    read's next base, not its anchor. Otherwise the read is not informative at any base of that
    allele. A read is never counted for an allele that starts with an indel with no aligned base
    before it, ends with a deletion with no aligned base after it, or holds a reference base other
    than A, C, G, or T, before or after left-alignment. Nor is it counted for an indel that would be
    left-aligned past the read's previous difference or reference skip, or past its first aligned
    base, so one read is never counted for two alleles at one base; it is then not informative from
    there to the end of the indel. A read with no stored qualities (QUAL `*`) has quality 255 at
    every base, as in htslib, so it passes every floor, and a read with no stored bases (SEQ `*`) is
    not counted at all.

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
                a `StreamingPileupBuilder` filters at 13 by default, as pysam's `pileup()` and
                `samtools mpileup` do.
            min_mapping_quality: the lowest mapping quality of a counted read.
            exclude_flags: reads with any of these SAM flags are not counted: by default,
                secondary, QC-fail, duplicate, and supplementary reads. Unmapped reads never are.

        Raises:
            ValueError: if a quality is not from 0 to 255, or `exclude_flags` from 0 to 65535.
        """
        self._reference: FastaFile = reference
        self._min_base_quality: int = min_base_quality
        self._min_mapping_quality: int = min_mapping_quality
        self._exclude_flags: int = exclude_flags
        self._tabulation: _native.Tabulation = _native.Tabulation(
            reference,
            min_base_quality=min_base_quality,
            min_mapping_quality=min_mapping_quality,
            exclude_flags=exclude_flags,
        )

    @property
    def reference(self) -> FastaFile:
        """The indexed reference the reads are aligned to."""
        return self._reference

    @property
    def min_base_quality(self) -> int:
        """The lowest base quality of an informative base."""
        return self._min_base_quality

    @property
    def min_mapping_quality(self) -> int:
        """The lowest mapping quality of a counted read."""
        return self._min_mapping_quality

    @property
    def exclude_flags(self) -> int:
        """The SAM flags of reads that are not counted."""
        return self._exclude_flags

    def accepts(self, record: AlignedSegment) -> bool:
        """Whether a read passes the flag and mapping-quality filters and has bases to count."""
        return self._tabulation.accepts(record)

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
        contigs = _by_header(alignments, self._reference, territory)
        return (
            base
            for contig, spans in contigs
            for base in self._tabulation.tabulate_contig(alignments, contig, spans)
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
        counted, dropped = self._tabulation.alleles(record)
        alleles = [
            Allele(pos=pos, ref=ref, alt=alt, start=start, end=end)
            for pos, ref, alt, start, end in counted
        ]
        return alleles, dropped


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
            `StreamingPileupBuilder` filters at 13 by default, as pysam's `pileup()` and
            `samtools mpileup` do.
        min_mapping_quality: the lowest mapping quality of a counted read.
        exclude_flags: reads with any of these SAM flags are not counted: by default, secondary,
            QC-fail, duplicate, and supplementary reads.

    Raises:
        ValueError: if a quality is not from 0 to 255, or `exclude_flags` from 0 to 65535, or a
            contig of the territory is not in the alignment header or the reference, or has
            another length in each, or a span runs past the end of its contig.
    """
    tabulator = Tabulator(
        reference,
        min_base_quality=min_base_quality,
        min_mapping_quality=min_mapping_quality,
        exclude_flags=exclude_flags,
    )
    return tabulator.tabulate(alignments, territory)
