from enum import StrEnum
from enum import auto
from typing import Final

from streampile import _native

DEFAULT_MIN_BASE_QUALITY: Final[int] = _native.DEFAULT_MIN_BASE_QUALITY
"""The default minimum base quality of a pileup, as in pysam's `pileup()` and `samtools mpileup`."""

DEFAULT_EXCLUDE_FLAGS: Final[int] = _native.DEFAULT_EXCLUDE_FLAGS
"""Secondary, QC-fail, duplicate, and supplementary reads, which are left out by default."""


class PileupReadType(StrEnum):
    """What a read holds at a pileup's position.

    A `skip` is a reference skip, the CIGAR `N` operator, not an `N` base, which is a `base`.
    """

    base = auto()
    deletion = auto()
    insertion = auto()
    skip = auto()


BASE = PileupReadType.base
DELETION = PileupReadType.deletion
INSERTION = PileupReadType.insertion
SKIP = PileupReadType.skip
