#!/usr/bin/env python3
"""Thin wrapper for users who prefer scripts/ entry points."""

from swarm_comsol_export.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
