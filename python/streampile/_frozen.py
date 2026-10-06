from collections.abc import Callable
from dataclasses import dataclass
from typing import TypeVar

from typing_extensions import dataclass_transform

T = TypeVar("T")


@dataclass_transform(frozen_default=True)
def frozen(*, kw_only: bool = False, slots: bool = False) -> Callable[[type[T]], type[T]]:
    """Make a dataclass that type checkers treat as frozen.

    At runtime it is a plain dataclass, hashed by its fields as a frozen one is, since a frozen
    dataclass sets every field through `object.__setattr__` and takes two to three times as long
    to build.
    """

    def decorate(cls: type[T]) -> type[T]:
        return dataclass(cls, kw_only=kw_only, slots=slots, unsafe_hash=True)

    return decorate
