"""Common type aliases and literals used throughout the recount3 API.

Attributes:
    CacheMode: Literal ``"enable" | "disable" | "update"``. Controls how
        :meth:`~recount3.resource.R3Resource.download` interacts with the
        on-disk cache:

        * ``"enable"``: use the existing cached file; download only if
          missing (default).
        * ``"disable"``: bypass the cache entirely and stream the file
          directly to the destination.
        * ``"update"``: force a fresh download even if a cached copy
          already exists, then cache the result.

    CompatibilityMode: Literal ``"family" | "feature"``. Controls validation
        in :meth:`~recount3.bundle.R3ResourceBundle.stack_count_matrices`:

        * ``"family"``: all count resources must belong to the same
          high-level family (gene/exon versus junctions).
        * ``"feature"``: stricter: all resources must additionally share
          an identical feature space (same genomic unit and annotation).

    StringOrIterable: ``str | Iterable[str]``. Most search functions accept
        either a single string or an iterable of strings for each parameter.
        When an iterable is passed, the function computes the Cartesian
        product across all parameters and returns one
        :class:`~recount3.resource.R3Resource` per combination.

    FieldSpec: ``StringOrIterable | Callable[[Any], bool] | None``. The
        filter predicate accepted by
        :meth:`~recount3.bundle.R3ResourceBundle.filter` and
        :func:`~recount3.search.match_spec`. Three forms are accepted:

        * ``None``: no filtering; every value passes.
        * A string or iterable of strings: exact membership test against the
          field value.
        * A callable: called with the field value; a truthy return keeps the
          resource.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from typing import Any, Literal, TypeAlias

CacheMode: TypeAlias = Literal["enable", "disable", "update"]
CompatibilityMode: TypeAlias = Literal["family", "feature"]

StringOrIterable: TypeAlias = str | Iterable[str]
FieldSpec: TypeAlias = StringOrIterable | Callable[[Any], bool] | None
