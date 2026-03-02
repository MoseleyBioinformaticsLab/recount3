#!/usr/bin/env python3
# pylint: disable=invalid-name
"""Module entry point.

Delegates to :mod:`recount3.cli`.
"""

from __future__ import annotations

from recount3.cli import main

if __name__ == "__main__":
    main()
