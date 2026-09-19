#!/usr/bin/env python3
"""Run the patched MCS, charge-renormalization, and parameter build for the pair."""

from __future__ import annotations

import contextlib
import json
import os
import shutil
import sys
import traceback
from pathlib import Path


HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"
PREP_ROOT = HERE.parents[1]
sys.path.insert(0, str(PREP_ROOT))

import msld_crn  # noqa: E402
import msld_mcs_rdecomp as mcs  # noqa: E402
import msld_prm  # noqa: E402


@contextlib.contextmanager
def working_directory(path: Path):
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def prepare(run_dir: Path) -> None:
    if run_dir.exists():
        raise RuntimeError(f"Refusing to overwrite existing run directory: {run_dir}")
    run_dir.mkdir(parents=True)
    stems = ["meso_erythritol", "threitol"]
    (run_dir / "mol_list.txt").write_text("\n".join(stems) + "\n")
    for stem in stems:
        for suffix in ("sdf", "mol2", "rtf", "prm", "str"):
            shutil.copyfile(INPUTS / f"{stem}.{suffix}", run_dir / f"{stem}.{suffix}")

    status = {"ligands": stems, "mcs": "not-run", "build": "not-run"}
    original_update = mcs.update_atomtypes
    mcs.update_atomtypes = lambda atomtypes: None
    try:
        with (run_dir / "run.log").open("w") as log:
            with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
                try:
                    with working_directory(run_dir):
                        mcs.MCSS_RDecomp("mol_list.txt", "MCS_for_MSLD.txt")
                    status["mcs"] = "success"
                    with working_directory(run_dir):
                        msld_crn.MsldCRN(
                            "MCS_for_MSLD.txt",
                            "build",
                            [[]],
                            [[]],
                            ChkQChange=True,
                            verbose=True,
                            debug=False,
                            ll=False,
                        )
                        msld_prm.MsldPRM(
                            "build", cgenff=True, verbose=True, debug=False
                        )
                    status["build"] = "success"
                except Exception as error:
                    status["error"] = f"{type(error).__name__}: {error}"
                    traceback.print_exc()
                    raise
    finally:
        mcs.update_atomtypes = original_update
        (run_dir / "status.json").write_text(json.dumps(status, indent=2) + "\n")


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: prepare_pair.py RUN_DIRECTORY")
    prepare(Path(sys.argv[1]).resolve())


if __name__ == "__main__":
    main()
