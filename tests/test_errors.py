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