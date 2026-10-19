from collections.abc import Iterable
from collections.abc import Iterator
from typing import TextIO

from pysam import AlignmentFile

from streampile import _native
from streampile._frozen import frozen
from streampile._pileup import DEFAULT_EXCLUDE_FLAGS


@frozen(slots=True)
class ContigBins:
    """The templates that start in each window of one contig.

    Attributes:
        contig: the name of the contig.
        window: the width of each window: the window asked for, or the length of the contig when
            it is shorter, as hmmcopy's `readCounter`, `gcCounter`, and `mapCounter` write them.
        templates: the templates that start in each window, from the start of the contig.
    """

    contig: str
    window: int
    templates: tuple[int, ...]


def bin_templates(
    alignments: AlignmentFile,
    window: int,
    *,
    contigs: Iterable[str] | None = None,
    min_mapping_quality: int = 0,
    exclude_flags: int = DEFAULT_EXCLUDE_FLAGS,
) -> Iterator[ContigBins]:
    """Iterate over contigs with the templates that start in each of their windows.

    A template is counted once, in the window where the first of its reads in coordinate order to
    pass the flag and mapping-quality filters starts. So the two reads of a pair count once, and a
    pair with one read filtered out still counts once, at the other. Mates are found by name, at
    the contig and position each read gives for its mate, its `RNEXT` and `PNEXT`, and a read with
    no name, or whose mate is unmapped or on a contig that is not binned, is counted alone.
    Unmapped, secondary, and supplementary records are never counted, so `exclude_flags` must hold
    the secondary and supplementary flags.

    Windows are 0-based and half-open from the start of each contig, so a contig has its length
    divided by the window, rounded up, of them, and the last is shorter when the window does not
    divide the length.

    Contigs come in the order of `contigs`, or of the alignment header, and are checked when this
    is called, before any read is read.

    Args:
        alignments: an indexed, coordinate-sorted BAM or CRAM.
        window: the width of each window.
        contigs: the contigs to bin, in order: by default, every contig of the alignment header.
        min_mapping_quality: the lowest mapping quality of a counted read.
        exclude_flags: reads with any of these SAM flags are not counted: by default, secondary,
            QC-fail, duplicate, and supplementary reads.

    Raises:
        ValueError: if a contig is not in the alignment header or is given twice, `window` is not
            from 1 to 4294967295, `min_mapping_quality` from 0 to 255, or `exclude_flags` from 0
            to 65535 with the secondary and supplementary flags.
    """
    lengths: dict[str, int] = {}
    for contig in alignments.references if contigs is None else contigs:
        if alignments.get_tid(contig) < 0:
            raise ValueError(f"Contig {contig} is not in the alignment header.")
        if contig in lengths:
            raise ValueError(f"Contig {contig} is given more than once.")
        lengths[contig] = alignments.get_reference_length(contig)
    bins = _native.TemplateBins(
        alignments,
        [(contig, alignments.get_tid(contig), length) for contig, length in lengths.items()],
        window,
        min_mapping_quality=min_mapping_quality,
        exclude_flags=exclude_flags,
    )
    return (
        ContigBins(contig=contig, window=min(window, lengths[contig]), templates=tuple(templates))
        for contig, templates in bins
    )


def write_wig(bins: Iterable[ContigBins], handle: TextIO) -> None:
    """Write bins as a fixedStep WIG, one block per contig, as hmmcopy's `readCounter` does.

    Args:
        bins: the bins of each contig, in order.
        handle: the text stream to write to.
    """
    for contig in bins:
        _ = handle.write(
            f"fixedStep chrom={contig.contig} start=1 step={contig.window} span={contig.window}\n"
        )
        handle.writelines(f"{templates}\n" for templates in contig.templates)
