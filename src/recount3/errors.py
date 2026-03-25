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
