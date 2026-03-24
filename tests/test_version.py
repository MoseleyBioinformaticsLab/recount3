from __future__ import annotations

import re

import recount3.version


def test_version_is_string() -> None:
    """__version__ must be a non-empty string."""
    assert isinstance(recount3.version.__version__, str)
    assert recount3.version.__version__


def test_version_matches_expected_value() -> None:
    """The pinned version value should match the package source."""
    assert recount3.version.__version__ == "1.0.0"


def test_version_is_semver_like() -> None:
    """__version__ should be in a simple MAJOR.MINOR.PATCH form."""
    pattern = re.compile(r"^[0-9]+\.[0-9]+\.[0-9]+$")
    assert pattern.fullmatch(recount3.version.__version__) is not None