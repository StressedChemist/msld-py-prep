#!/usr/bin/env python3
"""Run a solvated conventional CHARMM MSλD stereochemistry test.

The script reuses the minimized TIP3P box from ``solvated_lambda_test.py``,
builds a CHARMM PSF for those waters, runs continuous two-state theta/MSλD,
and measures the signed C3 volumes in every saved coordinate frame.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
from pathlib import Path

import mdtraj as md
import numpy as np


LIGAND_ATOMS = 22
EXPECTED_SIGNS = {"meso": -1, "threo": 1}


def read_pdb_atoms(filename: Path):
    atoms = []
    box_angstrom = None
    for line in filename.read_text().splitlines():
        if line.startswith("CRYST1"):
            box_angstrom = float(line[6:15])
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atoms.append(
            {
                "name": line[12:16].strip(),
                "resname": line[17:21].strip(),
                "resid": int(line[22:26]),
                "xyz": np.array(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                ),
                "element": line[76:78].strip(),
            }
        )
    if box_angstrom is None:
        raise ValueError(f"{filename} has no CRYST1 box record")
    if len(atoms) <= LIGAND_ATOMS:
        raise ValueError(f"{filename} contains no solvent")
    return atoms, box_angstrom


def pdb_record(serial, name, resname, resid, xyz, segid, element):
    return (
        f"ATOM  {serial:5d} {name:<4s} {resname:<4s}{resid:5d}    "
        f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
        f"  1.00  0.00      {segid:<4s}{element:>2s}\n"
    )


def write_charmm_coordinates(source: Path, output_dir: Path):
    atoms, box_angstrom = read_pdb_atoms(source)
    ligand_centre = np.mean([atom["xyz"] for atom in atoms[:LIGAND_ATOMS]], axis=0)
    for atom in atoms:
        atom["xyz"] -= ligand_centre

    ligand_lines = []
    for serial, atom in enumerate(atoms[:LIGAND_ATOMS], start=1):
        ligand_lines.append(
            pdb_record(
                serial,
                atom["name"],
                "LIG",
                1,
                atom["xyz"],
                "LIG",
                atom["element"],
            )
        )
    ligand_lines.extend(("TER\n", "END\n"))
    (output_dir / "ligand.pdb").write_text("".join(ligand_lines))

    waters = atoms[LIGAND_ATOMS:]
    if len(waters) % 3:
        raise ValueError("water atom count is not divisible by three")
    water_lines = []
    atom_names = ("OH2", "H1", "H2")
    for atom_index, atom in enumerate(waters):
        water_lines.append(
            pdb_record(
                atom_index + 1,
                atom_names[atom_index % 3],
                "TIP3",
                atom_index // 3 + 1,
                atom["xyz"],
                "WAT",
                atom["element"],
            )
        )
    water_lines.extend(("TER\n", "END\n"))
    (output_dir / "water.pdb").write_text("".join(water_lines))
    return box_angstrom, len(waters) // 3


def stage_charmm_inputs(run_dir: Path, output_dir: Path):
    prep_root = run_dir.parents[3]
    toppar = prep_root / "toppar"
    build = run_dir / "build"
    sources = {
        toppar / "top_all36_cgenff.rtf": "top_all36_cgenff.rtf",
        toppar / "par_all36_cgenff.prm": "par_all36_cgenff.prm",
        toppar / "toppar_water_ions.str": "toppar_water_ions.str",
        build / "full_ligand.prm": "full_ligand.prm",
        run_dir / "hybrid.psf": "hybrid.psf",
    }
    for source, destination in sources.items():
        shutil.copyfile(source, output_dir / destination)


def charmm_input(
    box_angstrom: float,
    steps: int,
    save_frequency: int,
    meso_bias: float,
    threo_bias: float,
    lambda_temperature: float,
    lambda_mass: float,
    lambda_friction: float,
):
    return f"""* Solvated meso-erythritol/threitol continuous MSLD test
*

bomblev -2

open read card unit 10 name top_all36_cgenff.rtf
read rtf card unit 10
open read card unit 20 name par_all36_cgenff.prm
read para card unit 20 flex
stream toppar_water_ions.str
read param flex append card name full_ligand.prm

read psf card name hybrid.psf
read coor pdb resid name ligand.pdb

read sequ pdb name water.pdb
generate WAT first none last none setup noangle nodihedral
read coor pdb resid name water.pdb

print coor sele .not. init end
write psf card name solvated.psf
* solvated dual-topology stereoisomer system
*
write coor pdb card name solvated_initial.pdb
* initial coordinates
*

crystal define cubic {box_angstrom:.6f} {box_angstrom:.6f} {box_angstrom:.6f} 90.0 90.0 90.0
crystal build cutoff 9.0 nope 0
image byres xcen 0 ycen 0 zcen 0 sele resn TIP3 end
image byres xcen 0 ycen 0 zcen 0 sele segid LIG end

define site1_sub1 select atom LIG 1 O015 .or. atom LIG 1 H016 .or. -
                         atom LIG 1 H017 .or. atom LIG 1 C018 end
define site1_sub2 select atom LIG 1 O019 .or. atom LIG 1 H020 .or. -
                         atom LIG 1 H021 .or. atom LIG 1 C022 end

block 3
   call 2 sele site1_sub1 show end
   call 3 sele site1_sub2 show end
   qldm theta
   lang temp {lambda_temperature:.6f}
   soft on
   ldin 1 1.0 0.0 {lambda_mass:.6f} 0.0 {lambda_friction:.6f}
   ldin 2 0.5 0.0 {lambda_mass:.6f} {meso_bias:.6f} {lambda_friction:.6f}
   ldin 3 0.5 0.0 {lambda_mass:.6f} {threo_bias:.6f} {lambda_friction:.6f}
   excl 2 3
   rmla bond thet impr
   msld 0 1 1 f2si
   msma
end

nbonds atom vatom vfswitch bycb cdie eps 1.0 -
   cutnb 9.0 ctofnb 8.0 ctonnb 7.0 cutim 9.0 -
   e14fac 1.0 wmin 1.5 imgfrq 20 inbfrq 20

energy
mini sd nstep 100 nprint 20
mini abnr nstep 400 nprint 50

shake fast para bonh tol 1.0e-8
scalar fbeta set 10.0 sele .not. hydrogen end

ldtitle
* meso/threo continuous MSλD lambda trajectory
*
open unit 21 write unform name trajectory.dcd
open unit 22 write form name restart.res
open unit 24 write unform name lambda.lmd

dynamics lang leap start timestep 0.001 nstep {steps} nprint 1000 iprfrq 1000 -
   iseed 34117 firstt 300.0 finalt 300.0 teminc 0.0 -
   ihtfrq 0 ieqfrq 0 ntrfrq 0 -
   iasors 1 iasvel 1 iscvel 0 ichecw 0 twindh 5.0 twindl -5.0 -
   iunread -1 iunwri 22 iuncrd 21 nsavc {save_frequency} -
   iunldm 24 nsavl {save_frequency} nsavv 0 iunvel -1 -
   eps 1.0 echeck 0

write coor pdb card name solvated_final.pdb
* final coordinates
*

stop
"""


def signed_volume(positions, names, centre_name, oxygen_name):
    centre = positions[names[centre_name]]
    return float(
        np.dot(
            np.cross(
                positions[names["C002"]] - centre,
                positions[names["C009"]] - centre,
            ),
            positions[names[oxygen_name]] - centre,
        )
    )


def parse_lambda_output(log_text: str):
    """Extract block-2/block-3 lambdas printed by CHARMM TRAJ LAMB."""
    rows = []
    pattern = re.compile(
        r"^\s*LAMBDA>\s+([-+0-9.EeDd]+)\s+"
        r"([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s*$"
    )
    for line in log_text.splitlines():
        match = pattern.match(line)
        if not match:
            continue
        values = [
            float(value.replace("D", "E").replace("d", "e"))
            for value in match.groups()
        ]
        rows.append(tuple(values))
    return rows


def lambda_motion_metrics(lambda_meso):
    dominant = (lambda_meso >= 0.5).astype(int)
    endpoint_sequence = []
    for value in lambda_meso:
        endpoint = 0 if value <= 0.1 else 1 if value >= 0.9 else None
        if endpoint is not None and (
            not endpoint_sequence or endpoint != endpoint_sequence[-1]
        ):
            endpoint_sequence.append(endpoint)
    return {
        "dominant_state_transitions": int(np.count_nonzero(np.diff(dominant))),
        "endpoint_alternations_lambda_0.1_0.9": max(0, len(endpoint_sequence) - 1),
        "meso_endpoint_frames_lambda_ge_0.9": int(
            np.count_nonzero(lambda_meso >= 0.9)
        ),
        "threo_endpoint_frames_lambda_le_0.1": int(
            np.count_nonzero(lambda_meso <= 0.1)
        ),
        "lambda_meso_range": [float(lambda_meso.min()), float(lambda_meso.max())],
    }


def analyze(
    output_dir: Path,
    save_frequency: int,
    minimum_abs_volume_nm3: float = 0.001,
):
    trajectory = md.load_dcd(
        str(output_dir / "trajectory.dcd"), top=str(output_dir / "solvated.psf")
    )
    ligand_atoms = [
        atom for atom in trajectory.topology.atoms if atom.residue.name == "LIG"
    ]
    names = {atom.name: atom.index for atom in ligand_atoms}
    if len(names) != len(ligand_atoms):
        raise RuntimeError("ligand atom names are not unique in the CHARMM topology")
    required_names = {"C002", "C009", "C018", "O015", "C022", "O019"}
    missing_names = sorted(required_names - names.keys())
    if missing_names:
        raise RuntimeError(
            "CHARMM topology lacks signed-volume atoms: " + ", ".join(missing_names)
        )
    meso = np.array(
        [signed_volume(frame, names, "C018", "O015") for frame in trajectory.xyz]
    )
    threo = np.array(
        [signed_volume(frame, names, "C022", "O019") for frame in trajectory.xyz]
    )
    lambdas = parse_lambda_output((output_dir / "lambda.log").read_text())
    if len(lambdas) != trajectory.n_frames:
        raise RuntimeError(
            f"coordinate/lambda frame mismatch: {trajectory.n_frames} versus "
            f"{len(lambdas)}"
        )
    if trajectory.n_frames < 2:
        raise RuntimeError("at least two synchronized frames are required")

    lambda_values = np.array([[row[1], row[2]] for row in lambdas])
    if np.any(lambda_values < -0.0001) or np.any(lambda_values > 1.0001):
        raise RuntimeError("CHARMM lambda values fall outside [0, 1]")
    maximum_lambda_sum_deviation = float(
        np.max(np.abs(np.sum(lambda_values, axis=1) - 1.0))
    )
    if maximum_lambda_sum_deviation > 0.001:
        raise RuntimeError("the two substituent lambdas do not sum to one")

    rows = []
    for frame_index in range(trajectory.n_frames):
        lambda_meso = lambdas[frame_index][1]
        lambda_threo = lambdas[frame_index][2]
        rows.append(
            {
                "frame": frame_index + 1,
                "step": (frame_index + 1) * save_frequency,
                "time_ps": lambdas[frame_index][0],
                "lambda_meso": lambda_meso,
                "lambda_threo": lambda_threo,
                "meso_signed_volume_nm3": meso[frame_index],
                "threo_signed_volume_nm3": threo[frame_index],
            }
        )
    with (output_dir / "lambda_stereo.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    lambda_meso = np.array([row[1] for row in lambdas])
    midpoint = len(lambda_meso) // 2
    full_motion = lambda_motion_metrics(lambda_meso)
    first_half_motion = lambda_motion_metrics(lambda_meso[:midpoint])
    second_half_motion = lambda_motion_metrics(lambda_meso[midpoint:])
    motion_persisted = all(
        metrics["dominant_state_transitions"] > 0
        and metrics["endpoint_alternations_lambda_0.1_0.9"] > 0
        for metrics in (first_half_motion, second_half_motion)
    )
    separated = bool(
        np.all(meso < -minimum_abs_volume_nm3)
        and np.all(threo > minimum_abs_volume_nm3)
    )
    summary = {
        "coordinate_frames": trajectory.n_frames,
        "lambda_frames": len(lambdas),
        "maximum_lambda_sum_deviation": maximum_lambda_sum_deviation,
        "lambda_motion": full_motion,
        "lambda_motion_first_half": first_half_motion,
        "lambda_motion_second_half": second_half_motion,
        "lambda_motion_persisted": motion_persisted,
        "meso_signed_volume_nm3": [float(meso.min()), float(meso.max())],
        "threo_signed_volume_nm3": [float(threo.min()), float(threo.max())],
        "meso_sign_flips": int(np.count_nonzero(np.sign(meso) != EXPECTED_SIGNS["meso"])),
        "threo_sign_flips": int(np.count_nonzero(np.sign(threo) != EXPECTED_SIGNS["threo"])),
        "minimum_required_abs_signed_volume_nm3": minimum_abs_volume_nm3,
        "separated_throughout": separated,
        "working_construction": bool(separated and motion_persisted),
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def run_lambda_analysis(charmm: Path, output_dir: Path, environment):
    analysis_input = """* Print continuous MSLD lambda trajectory
*
bomblev -2
open unit 24 read unform name lambda.lmd
traj lamb print first 24 nunit 1 nosub
close unit 24
stop
"""
    (output_dir / "lambda_analysis.inp").write_text(analysis_input)
    with (output_dir / "lambda_analysis.inp").open("rb") as input_handle, (
        output_dir / "lambda.log"
    ).open("wb") as output_handle:
        subprocess.run(
            [str(charmm)],
            stdin=input_handle,
            stdout=output_handle,
            stderr=subprocess.STDOUT,
            cwd=output_dir,
            env=environment,
            check=False,
        )
    rows = parse_lambda_output((output_dir / "lambda.log").read_text(errors="replace"))
    if not rows:
        raise RuntimeError("CHARMM lambda analysis produced no readable frames")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--output-name", default="charmm_msld")
    parser.add_argument("--charmm", type=Path, default=shutil.which("charmm"))
    parser.add_argument("--steps", type=int, default=50000)
    parser.add_argument("--save-frequency", type=int, default=10)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--meso-bias", type=float, default=0.0)
    parser.add_argument("--threo-bias", type=float, default=0.0)
    parser.add_argument("--lambda-temperature", type=float, default=300.0)
    parser.add_argument("--lambda-mass", type=float, default=5.0)
    parser.add_argument("--lambda-friction", type=float, default=5.0)
    parser.add_argument("--minimum-abs-volume-nm3", type=float, default=0.001)
    parser.add_argument("--require-pass", action="store_true")
    args = parser.parse_args()
    if args.charmm is None:
        raise RuntimeError("CHARMM is not on PATH; load it or pass --charmm")

    run_dir = args.run_dir.resolve()
    output_dir = run_dir / args.output_name
    if output_dir.exists():
        raise RuntimeError(f"Refusing to overwrite existing directory: {output_dir}")
    output_dir.mkdir(parents=True)

    source = run_dir / "solvated_lambda" / "solvated_minimized.pdb"
    box_angstrom, water_count = write_charmm_coordinates(source, output_dir)
    stage_charmm_inputs(run_dir, output_dir)
    (output_dir / "run.inp").write_text(
        charmm_input(
            box_angstrom,
            args.steps,
            args.save_frequency,
            args.meso_bias,
            args.threo_bias,
            args.lambda_temperature,
            args.lambda_mass,
            args.lambda_friction,
        )
    )
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = str(args.threads)
    with (output_dir / "run.inp").open("rb") as input_handle, (
        output_dir / "charmm.log"
    ).open("wb") as output_handle:
        result = subprocess.run(
            [str(Path(args.charmm).resolve())],
            stdin=input_handle,
            stdout=output_handle,
            stderr=subprocess.STDOUT,
            cwd=output_dir,
            env=environment,
            check=False,
        )
    if result.returncode:
        raise RuntimeError(f"CHARMM exited with status {result.returncode}; inspect charmm.log")
    log_text = (output_dir / "charmm.log").read_text(errors="replace")
    if "ABNORMAL TERMINATION" in log_text:
        raise RuntimeError("CHARMM reported abnormal termination; inspect charmm.log")
    version_line = next(
        (line.strip() for line in log_text.splitlines() if "(CHARMM)" in line),
        "unknown",
    )
    commit_line = next(
        (line.strip() for line in log_text.splitlines() if "Git commit ID" in line),
        "unknown",
    )

    run_lambda_analysis(Path(args.charmm).resolve(), output_dir, environment)
    summary = analyze(
        output_dir,
        args.save_frequency,
        minimum_abs_volume_nm3=args.minimum_abs_volume_nm3,
    )
    summary["conditions"] = {
        "engine": "CHARMM",
        "engine_version": version_line,
        "engine_commit": commit_line,
        "lambda_method": "continuous theta/MSLD F2SI",
        "solvent": "TIP3P",
        "water_molecules": water_count,
        "box_angstrom": box_angstrom,
        "steps": args.steps,
        "timestep_fs": 1.0,
        "save_frequency": args.save_frequency,
        "temperature_K": 300.0,
        "ensemble": "NVT",
        "nonbonded_method": "periodic cutoff with force switching",
        "cutoff_angstrom": 8.0,
        "meso_fixed_bias_kcal_mol": args.meso_bias,
        "threo_fixed_bias_kcal_mol": args.threo_bias,
        "lambda_temperature_K": args.lambda_temperature,
        "lambda_mass": args.lambda_mass,
        "lambda_friction": args.lambda_friction,
        "random_seed": 34117,
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    if args.require_pass and not summary["working_construction"]:
        raise RuntimeError("MSLD construction validation did not pass")


if __name__ == "__main__":
    main()
