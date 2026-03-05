"""Common type aliases and literals."""

from __future__ import annotations

from collections.abc import Callable, Iterable
from typing import Any, Literal, TypeAlias

CacheMode: TypeAlias = Literal["enable", "disable", "update"]
CompatibilityMode: TypeAlias = Literal["family", "feature"]

StringOrIterable: TypeAlias = str | Iterable[str]
FieldSpec: TypeAlias = StringOrIterable | Callable[[Any], bool] | None
