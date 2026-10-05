from os import PathLike
from typing import final

from typing_extensions import Self

NativeEntry = tuple[
    str,
    int,
    str,
    int | None,
    int | None,
    int | None,
    int,
    str | None,
    int | None,
    str | None,
    bytes | None,
]
"""An entry: read name, flag, kind, query position, query position or next, insertion offset,
insertion length, base, quality, inserted bases, and inserted qualities."""

NativeColumn = tuple[list[NativeEntry], int, int, str, bytes]
"""A column: its entries, unfiltered depth, filtered depth, bases, and qualities."""

def sweep(
    path: str | PathLike[str],
    spans: list[tuple[str, int, int]],
    *,
    min_mapping_quality: int = 0,
    exclude_flags: int = 0xF00,
    min_base_quality: int = 13,
    proper_pairs_only: bool = False,
    without_overlaps: bool = False,
) -> list[NativeColumn]:
    """Sweep spans of a BAM with the Rust builder and return every column."""

@final
class ColumnCounts:
    """Count the bases, deletions, and insertions of every column of spans of a BAM."""

    def __new__(
        cls,
        path: str | PathLike[str],
        spans: list[tuple[str, int, int]],
        *,
        min_mapping_quality: int,
        exclude_flags: int,
        quality_floor: int,
    ) -> Self: ...
    def __iter__(self) -> Self: ...
    def __next__(self) -> dict[str, int]: ...
