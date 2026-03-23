#!/usr/bin/env python3
# pylint: disable=invalid-name
"""Module entry point.

Delegates to :mod:`recount3.cli`. Invoke via ``python -m recount3 --help``.
"""

from __future__ import annotations

from recount3.cli import main

if __name__ == "__main__":
    main()
