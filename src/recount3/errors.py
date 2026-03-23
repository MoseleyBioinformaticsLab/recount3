"""Domain-specific exceptions for recount3.

All exceptions inherit from :class:`Recount3Error`, so callers can catch
the base class to handle any recount3-specific failure, or catch a subclass
to handle a specific failure mode:

* :class:`Recount3Error`: base class for all recount3 exceptions.
* :class:`ConfigurationError`: invalid or missing configuration (bad env
  var values, inaccessible cache directory, unsupported option combinations).
* :class:`DownloadError`: a network or I/O failure occurred while
  downloading a resource.
* :class:`LoadError`: a resource was downloaded but could not be parsed
  (empty file, unexpected format, shape mismatch, missing columns).
* :class:`CompatibilityError`: resources that are incompatible with each
  other were combined in an operation such as
  :meth:`~recount3.bundle.R3ResourceBundle.stack_count_matrices`.

Example:
    Catch all recount3 errors with the base class::

        from recount3.errors import Recount3Error, DownloadError

        try:
            res.download(path="/data")
        except DownloadError as exc:
            print(f"Network failure: {exc}")
        except Recount3Error as exc:
            print(f"Unexpected recount3 error: {exc}")
"""

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
    """Raised when resources are incompatible for combined operations."""
