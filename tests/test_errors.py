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

from __future__ import annotations

import recount3.errors

import pytest


def _all_error_types() -> tuple[type[BaseException], ...]:
    """Return all public recount3 exception types in stable order."""
    return (
        recount3.errors.Recount3Error,
        recount3.errors.ConfigurationError,
        recount3.errors.DownloadError,
        recount3.errors.LoadError,
        recount3.errors.CompatibilityError,
    )


def test_error_types_are_exceptions() -> None:
    """All exported error types should be subclasses of Exception."""
    for error_type in _all_error_types():
        assert issubclass(error_type, Exception)


def test_error_hierarchy_is_correct() -> None:
    """All domain errors should derive from the package base error."""
    base_type = recount3.errors.Recount3Error
    leaf_types = (
        recount3.errors.ConfigurationError,
        recount3.errors.DownloadError,
        recount3.errors.LoadError,
        recount3.errors.CompatibilityError,
    )

    assert issubclass(base_type, Exception)
    for leaf_type in leaf_types:
        assert issubclass(leaf_type, base_type)
        assert leaf_type is not base_type


@pytest.mark.parametrize(
    "error_type",
    (
        recount3.errors.Recount3Error,
        recount3.errors.ConfigurationError,
        recount3.errors.DownloadError,
        recount3.errors.LoadError,
        recount3.errors.CompatibilityError,
    ),
)
def test_error_is_catchable_as_base(error_type: type[BaseException]) -> None:
    """Each specific error should be catchable via Recount3Error."""
    message = "unit-test-message"
    try:
        raise error_type(message)
    except recount3.errors.Recount3Error as caught:
        assert str(caught) == message


@pytest.mark.parametrize(
    "error_type",
    (
        recount3.errors.Recount3Error,
        recount3.errors.ConfigurationError,
        recount3.errors.DownloadError,
        recount3.errors.LoadError,
        recount3.errors.CompatibilityError,
    ),
)
def test_error_message_is_visible_in_repr(
    error_type: type[BaseException],
) -> None:
    """repr(exception) should include the provided message for debugging."""
    message = "debuggable"
    exc = error_type(message)
    assert message in repr(exc)