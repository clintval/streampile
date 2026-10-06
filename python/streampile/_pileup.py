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


class AgreementStrategy(StrEnum):
    """How the quality of a template is made from two of its reads holding the same base.

    These are the strategies of fgbio's `CallOverlappingConsensusBases`, as fgumi implements them.

    Attributes:
        consensus: the sum of the two qualities, at most 93.
        max_qual: the higher of the two qualities.
        pass_through: each read keeps its quality, so the template has the higher of the two.
    """

    consensus = auto()
    max_qual = auto()
    pass_through = auto()


class DisagreementStrategy(StrEnum):
    """How a template's base and quality are made from two of its reads holding different bases.

    These are the strategies of fgbio's `CallOverlappingConsensusBases`, as fgumi implements them.

    Attributes:
        consensus: the base of the higher quality, at the higher quality less the lower and at
            least 2, or an `N` at quality 2 when the qualities are equal.
        mask_both: an `N` at quality 2.
        mask_lower_qual: the base of the higher quality at its own quality, the other masked, or an
            `N` at quality 2 when the qualities are equal.
    """

    consensus = auto()
    mask_both = auto()
    mask_lower_qual = auto()


BASE = PileupReadType.base
DELETION = PileupReadType.deletion
INSERTION = PileupReadType.insertion
SKIP = PileupReadType.skip
