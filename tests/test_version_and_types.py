"""Import checks for version string and public enums."""

from __future__ import annotations

import re

from recount3.version import __version__
from recount3.types import CacheMode, CompatibilityMode


def test_version_semver_like() -> None:
    """`__version__` should look like semantic versioning."""
    assert isinstance(__version__, str)
    assert re.match(r"^\d+\.\d+\.\d+.*$", __version__) is not None


def test_enum_members_accept_cli_strings() -> None:
    """Enum-like type aliases should accept the CLI string values."""
    # Cache mode values accepted by the CLI.
    for s in ("enable", "disable", "update"):
        assert s in CacheMode.__args__  # typing.Literal tuple
    for s in ("family", "feature"):
        assert s in CompatibilityMode.__args__
