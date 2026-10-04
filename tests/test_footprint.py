import pytest

from streampile._footprint import DELETED_AT_END
from streampile._footprint import SKIPPED
from streampile._footprint import Footprint
from streampile._footprint import is_placed

from .records import record
from .records import unmapped


def next_after(offset: int) -> int:
    return -offset - 3


@pytest.mark.parametrize(
    "cigar,bases,offsets,insertions,first,end",
    [
        ("4M", "ACGT", [0, 1, 2, 3], {}, 10, 14),
        ("2S4M", "TTACGT", [2, 3, 4, 5], {}, 10, 14),
        ("3H4M2H", "ACGT", [0, 1, 2, 3], {}, 10, 14),
        ("2M2D2M", "ACGT", [0, 1, next_after(2), next_after(2), 2, 3], {}, 10, 16),
        ("2M2N2M", "ACGT", [0, 1, SKIPPED, SKIPPED, 2, 3], {}, 10, 16),
        ("2M1I2M", "ACGTA", [0, 1, 3, 4], {11: (2, 1)}, 10, 14),
        ("1I3M", "ACGT", [1, 2, 3], {9: (0, 1)}, 9, 13),
        ("2S1I3M", "TTACGT", [3, 4, 5], {9: (2, 1)}, 9, 13),
        ("3M1I", "ACGT", [0, 1, 2], {12: (3, 1)}, 10, 13),
        ("2M1D1I2M", "ACGTA", [0, 1, next_after(2), 3, 4], {12: (2, 1)}, 10, 15),
        ("2M1I1I2M", "ACGTAC", [0, 1, 4, 5], {11: (2, 2)}, 10, 14),
        ("3M2D", "ACG", [0, 1, 2, DELETED_AT_END, DELETED_AT_END], {}, 10, 15),
        ("3M2D1S", "ACGT", [0, 1, 2, next_after(3), next_after(3)], {}, 10, 15),
        ("1=1X2M", "ACGT", [0, 1, 2, 3], {}, 10, 14),
    ],
)
def test_footprint_walks_the_cigar_once(
    cigar: str,
    bases: str,
    offsets: list[int],
    insertions: dict[int, tuple[int, int]],
    first: int,
    end: int,
) -> None:
    footprint = Footprint(record("r", 10, cigar, bases))
    assert list(footprint.offsets) == offsets
    assert footprint.insertions == insertions
    assert footprint.first == first
    assert footprint.start == 10
    assert footprint.end == end


@pytest.mark.parametrize(
    "cigar,bases,placed",
    [
        ("4M", "ACGT", True),
        ("1D3M", "ACG", True),
        ("2S2N2M", "ACGT", True),
        ("1=1X2S", "ACGT", True),
        ("4D", "*", True),
        ("4S", "ACGT", False),
        ("4I", "ACGT", False),
        ("2S2I", "ACGT", False),
        ("4H4S", "ACGT", False),
        ("1S2I1S", "ACGT", False),
    ],
)
def test_is_placed_needs_a_reference_consuming_operator(
    cigar: str, bases: str, placed: bool
) -> None:
    assert is_placed(record("r", 10, cigar, bases)) is placed


def test_is_placed_needs_a_mapped_read() -> None:
    assert not is_placed(unmapped("u"))
    assert not is_placed(record("r", 10, "4M", "ACGT", flag=4))
