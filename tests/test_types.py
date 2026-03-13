from __future__ import annotations

import collections.abc
import types
import typing

import recount3.types


def test_cache_mode_is_literal_with_expected_values() -> None:
    """CacheMode must be a Literal with the three supported cache modes."""
    origin = typing.get_origin(recount3.types.CacheMode)
    args = typing.get_args(recount3.types.CacheMode)

    assert origin is typing.Literal
    assert set(args) == {"enable", "disable", "update"}


def test_compatibility_mode_is_literal_with_expected_values() -> None:
    """CompatibilityMode must be a Literal describing supported combine modes."""
    origin = typing.get_origin(recount3.types.CompatibilityMode)
    args = typing.get_args(recount3.types.CompatibilityMode)

    assert origin is typing.Literal
    assert set(args) == {"family", "feature"}


def test_string_or_iterable_is_union_of_str_and_iterable_str() -> None:
    """StringOrIterable must accept str or an iterable of str."""
    union_type = recount3.types.StringOrIterable
    origin = typing.get_origin(union_type)
    args = typing.get_args(union_type)

    assert origin is types.UnionType
    assert set(args) == {str, collections.abc.Iterable[str]}


def test_field_spec_contains_expected_union_members() -> None:
    """FieldSpec must accept None, a string/iterable, or a predicate callable."""
    field_spec = recount3.types.FieldSpec
    origin = typing.get_origin(field_spec)
    args = typing.get_args(field_spec)

    expected_callable = collections.abc.Callable[[typing.Any], bool]
    expected_args = {
        str,
        collections.abc.Iterable[str],
        expected_callable,
        type(None),
    }

    assert origin is types.UnionType
    assert set(args) == expected_args

    callable_origin = typing.get_origin(expected_callable)
    callable_args = typing.get_args(expected_callable)
    assert callable_origin is collections.abc.Callable
    assert len(callable_args) == 2
    assert callable_args[1] is bool