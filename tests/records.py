from array import array
from collections.abc import Iterable
from pathlib import Path

import pysam
from pysam import AlignedSegment
from pysam import AlignmentFile
from pysam import AlignmentHeader

DATA = Path(__file__).parent / "data"

HEADER = AlignmentHeader.from_text(
    "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000\n"
)


def record(
    name: str,
    start: int,
    cigar: str,
    bases: str,
    *,
    contig: str = "chr1",
    flag: int = 0,
    mapq: int = 60,
    quals: list[int] | None = None,
    header: AlignmentHeader = HEADER,
) -> AlignedSegment:
    """A mapped read built field by field."""
    read = AlignedSegment(header)
    read.query_name = name
    read.flag = flag
    read.reference_name = contig
    read.reference_start = start
    read.mapping_quality = mapq
    read.cigarstring = cigar
    read.query_sequence = bases
    read.query_qualities = array("B", [40] * len(bases) if quals is None else quals)
    return read


def unmapped(name: str, bases: str = "ACGT", header: AlignmentHeader = HEADER) -> AlignedSegment:
    """A read with no position at all."""
    read = AlignedSegment(header)
    read.query_name = name
    read.flag = 4
    read.reference_id = -1
    read.reference_start = -1
    read.query_sequence = bases
    read.query_qualities = array("B", [40] * len(bases))
    return read


def write_bam(
    path: Path, reads: Iterable[AlignedSegment], header: AlignmentHeader = HEADER
) -> Path:
    """Write reads, already in coordinate order, to an indexed BAM."""
    with AlignmentFile(str(path), "wb", header=header) as sink:
        for read in reads:
            sink.write(read)
    pysam.index(str(path))
    return path
