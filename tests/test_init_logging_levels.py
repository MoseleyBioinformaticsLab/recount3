"""Unit tests for CLI logging level initialization.

Test the convenience wrapper `_init_logging` by asserting that messages at
DEBUG/INFO/WARNING make it through based on `--quiet`/`--verbose`.
"""

from __future__ import annotations

import logging

import recount3.cli as cli


def _args(**over):
    base = dict(quiet=False, verbose=False)
    base.update(over)
    return type("Args", (), base)()


def test_logging_verbose_enables_debug(caplog) -> None:
    """`--verbose` should enable DEBUG visibility."""
    caplog.set_level(logging.CRITICAL)  # ensure overridden by init
    cli._init_logging(_args(verbose=True))  # type: ignore
    logger = logging.getLogger("recount3.test")
    logger.debug("debug-visible")
    assert "debug-visible" in caplog.text


def test_logging_quiet_hides_info(caplog) -> None:
    """`--quiet` should suppress INFO messages."""
    cli._init_logging(_args(quiet=True))  # type: ignore
    logger = logging.getLogger("recount3.test")
    logger.info("info-hidden")
    assert "info-hidden" not in caplog.text
    logger.warning("warn-shown")
    assert "warn-shown" in caplog.text
