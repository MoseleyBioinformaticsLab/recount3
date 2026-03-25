# Copyright (c) 2026, Alexander A. Alsalihi, Robert M. Flight,
# Hunter N.B. Moseley. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# * All advertising materials mentioning features or use of this software must
#   display the following acknowledgement: This product includes software
#   developed by the copyright holder.
# * Neither the name of the copyright holder nor the names of its contributors
#   may be used to endorse or promote products derived from this software without
#   specific prior written permission.
# * If the source code is used in a published work, then proper citation of the
#   source code must be included with the published work.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS AS IS AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
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
