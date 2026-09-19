#!/usr/bin/env python3
"""Run stock and chirality-aware msld-py-prep decompositions."""

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
RUNS = HERE / "runs"
PREP_ROOT = HERE.parents[1]
sys.path.insert(0, str(PREP_ROOT))

import msld_crn  # noqa: E402
import msld_mcs_rdecomp as mcs  # noqa: E402


SCENARIOS = {
    "pair": "mol_list_pair.txt",
    "corrected_original_series": "mol_list_corrected_original_series.txt",
    "five_state_series": "mol_list_series.txt",
}


@contextlib.contextmanager
def working_directory(path: Path):
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(previous)


def copy_inputs(run_dir: Path, list_name: str) -> list[str]:
    stems = [line.strip() for line in (INPUTS / list_name).read_text().splitlines() if line.strip()]
    (run_dir / "mol_list.txt").write_text("\n".join(stems) + "\n")
    for stem in stems:
        for suffix in ("sdf", "str", "rtf", "mol2"):
            shutil.copyfile(INPUTS / f"{stem}.{suffix}", run_dir / f"{stem}.{suffix}")
    return stems


def site_count(mcs_path: Path) -> int:
    for line in mcs_path.read_text().splitlines():
        if line.startswith("NSUBS"):
            return len(line.split()) - 1
    raise RuntimeError(f"No NSUBS record in {mcs_path}")


def run_case(scenario: str, list_name: str, chirality_aware: bool) -> None:
    mode = "chiral" if chirality_aware else "stock"
    run_dir = RUNS / f"{mode}_{scenario}"
    run_dir.mkdir(parents=True, exist_ok=False)
    stems = copy_inputs(run_dir, list_name)
    status = {
        "scenario": scenario,
        "mode": mode,
        "ligands": stems,
        "mcs": "not-run",
        "charge_renormalization": "not-run",
    }

    original_update = mcs.update_atomtypes
    original_find_mcs = mcs.rdFMCS.FindMCS
    mcs.update_atomtypes = lambda atomtypes: None

    if chirality_aware:
        def find_mcs_with_chirality(*args, **kwargs):
            kwargs["matchChiralTag"] = True
            return original_find_mcs(*args, **kwargs)

        mcs.rdFMCS.FindMCS = find_mcs_with_chirality

    with (run_dir / "run.log").open("w") as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        try:
            with working_directory(run_dir):
                mcs.MCSS_RDecomp("mol_list.txt", "MCS_for_MSLD.txt")
            status["mcs"] = "success"
        except Exception as error:  # Preserve the exact stock failure for inspection.
            status["mcs"] = f"{type(error).__name__}: {error}"
            traceback.print_exc()
        finally:
            mcs.update_atomtypes = original_update
            mcs.rdFMCS.FindMCS = original_find_mcs

        mcs_path = run_dir / "MCS_for_MSLD.txt"
        if mcs_path.exists():
            try:
                nsites = site_count(mcs_path)
                with working_directory(run_dir):
                    msld_crn.MsldCRN(
                        "MCS_for_MSLD.txt",
                        "build",
                        [[] for _ in range(nsites)],
                        [[] for _ in range(nsites)],
                        ChkQChange=True,
                        verbose=True,
                        debug=False,
                        ll=False,
                    )
                status["charge_renormalization"] = "success"
            except Exception as error:
                status["charge_renormalization"] = f"{type(error).__name__}: {error}"
                traceback.print_exc()

    (run_dir / "status.json").write_text(json.dumps(status, indent=2) + "\n")


def main() -> None:
    RUNS.mkdir(parents=True, exist_ok=True)
    existing = [RUNS / f"{mode}_{scenario}" for scenario in SCENARIOS for mode in ("stock", "chiral")]
    collisions = [path for path in existing if path.exists()]
    if collisions:
        joined = "\n".join(str(path) for path in collisions)
        raise RuntimeError(f"Refusing to overwrite existing analysis directories:\n{joined}")

    for scenario, list_name in SCENARIOS.items():
        run_case(scenario, list_name, chirality_aware=False)
        run_case(scenario, list_name, chirality_aware=True)


if __name__ == "__main__":
    main()
