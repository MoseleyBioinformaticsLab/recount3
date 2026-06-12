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
#   may be used to endorse or promote products derived from this software
#   without specific prior written permission.
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
"""Platform-conditional skip markers for recount3's optional dependencies.

recount3's core is pure Python and runs everywhere, but two of its optional
features are backed by compiled (C / C++) dependencies whose binary-wheel
availability differs by operating system. Continuous integration therefore
installs a different subset of extras on each OS (Linux: full; macOS: no bigwig;
Windows: core only), and the tests that exercise those features must be skipped
wherever their dependency cannot be installed, otherwise the suite would fail
on a platform where the feature is simply not available. This module centralizes
that skipping logic. Each optional stack gets one marker.
"""

from __future__ import annotations

import importlib.util

import pytest


def _all_installed(*module_names: str) -> bool:
    """Return True only if every named module is importable on this platform."""
    return all(
        importlib.util.find_spec(name) is not None for name in module_names
    )


HAS_PYBIGWIG = _all_installed("pyBigWig")
HAS_BIOCPY = _all_installed("summarizedexperiment", "genomicranges", "iranges")

_OPTIONAL_DEPENDENCY_MARKERS: dict[str, tuple[bool, str]] = {
    "requires_pybigwig": (
        HAS_PYBIGWIG,
        "pyBigWig is not installed; it has no Windows support and is omitted on "
        "macOS in CI (installed on Linux only).",
    ),
    "requires_biocpy": (
        HAS_BIOCPY,
        "BiocPy ranged stack (summarizedexperiment/genomicranges/iranges) is not "
        "installed; iranges has no Windows wheel and fails to build there "
        "(installed on Linux and macOS).",
    ),
}


def pytest_configure(config: pytest.Config) -> None:
    """Register the optional-dependency markers so ``--strict-markers`` passes."""
    config.addinivalue_line(
        "markers",
        "requires_pybigwig: test needs the optional pyBigWig dependency "
        "(unavailable on Windows; omitted on macOS in CI).",
    )
    config.addinivalue_line(
        "markers",
        "requires_biocpy: test needs the BiocPy ranged stack "
        "(summarizedexperiment/genomicranges/iranges), unavailable on Windows.",
    )


def pytest_collection_modifyitems(items: list[pytest.Item]) -> None:
    """Skip tests whose required optional dependency is absent on this platform."""
    skip_marks = {
        name: pytest.mark.skip(reason=reason)
        for name, (available, reason) in _OPTIONAL_DEPENDENCY_MARKERS.items()
        if not available
    }
    if not skip_marks:
        return
    for item in items:
        for name, skip_mark in skip_marks.items():
            if name in item.keywords:
                item.add_marker(skip_mark)
