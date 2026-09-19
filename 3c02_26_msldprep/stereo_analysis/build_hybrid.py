#!/usr/bin/env python3
"""Execute build_hybrid.inp through pyCHARMM in a prepared pair directory."""

from __future__ import annotations

import contextlib
import os
import sys
from pathlib import Path

from pycharmm import lingo


HERE = Path(__file__).resolve().parent


@contextlib.contextmanager
def working_directory(path: Path):
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: build_hybrid.py RUN_DIRECTORY")
    run_dir = Path(sys.argv[1]).resolve()
    lines = (HERE / "build_hybrid.inp").read_text().splitlines()
    # Top-level CHARMM title records and STOP are file-driver syntax rather
    # than commands accepted by pyCHARMM's lingo interface.
    commands = "\n".join(lines[3:-1])
    with working_directory(run_dir):
        lingo.charmm_script(commands)


if __name__ == "__main__":
    main()
