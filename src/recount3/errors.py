"""Domain-specific exceptions for recount3."""

from __future__ import annotations


class Recount3Error(Exception):
    """Base class for recount3-related errors."""


class ConfigurationError(Recount3Error):
    """Raised when configuration is invalid or inconsistent."""


class DownloadError(Recount3Error):
    """Raised when a resource fails to download."""


class LoadError(Recount3Error):
    """Raised when a resource fails to load or parse."""


class CompatibilityError(Recount3Error):
    """Raised when incompatible resources are selected for a combined operation."""
