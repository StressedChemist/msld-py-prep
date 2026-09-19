#!/usr/bin/env python3
"""Run fixed-lambda endpoint tests and monitor C3 configuration.

The inactive substituent has its charge, Lennard-Jones epsilon, and all
associated exceptions set to zero.  Bonded terms (including the signed
improper) remain active, matching the LaDyBUGS treatment of bonded forces.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from openmm import LangevinMiddleIntegrator, NonbondedForce, Platform, unit
from openmm.app import CharmmParameterSet, CharmmPsfFile, HBonds, NoCutoff, PDBFile, Simulation


ENDPOINTS = {
    "meso": {
        "lambda": [1.0, 0.0],
        "fragment": ["O015", "H016", "H017", "C018"],
        "inactive": ["O019", "H020", "H021", "C022"],
        "c3": "C018",
        "o3": "O015",
        "improper": ["C018", "C009", "C002", "O015"],
        "expected_volume_sign": -1,
    },
    "threo": {
        "lambda": [0.0, 1.0],
        "fragment": ["O019", "H020", "H021", "C022"],
        "inactive": ["O015", "H016", "H017", "C018"],
        "c3": "C022",
        "o3": "O019",
        "improper": ["C022", "C002", "C009", "O019"],
        "expected_volume_sign": 1,
    },
}


def torsion_degrees(points: np.ndarray) -> float:
    p0, p1, p2, p3 = points
    b0 = -(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return math.degrees(math.atan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))


def signed_volume(positions: np.ndarray, names: dict[str, int], endpoint: dict) -> float:
    centre = positions[names[endpoint["c3"]]]
    return float(
        np.dot(
            np.cross(
                positions[names["C002"]] - centre,
                positions[names["C009"]] - centre,
            ),
            positions[names[endpoint["o3"]]] - centre,
        )
    )


def deactivate_fragment(system, inactive_indices: set[int]) -> None:
    nonbonded = next(force for force in system.getForces() if isinstance(force, NonbondedForce))
    for atom_index in inactive_indices:
        charge, sigma, _ = nonbonded.getParticleParameters(atom_index)
        nonbonded.setParticleParameters(atom_index, 0.0 * charge, sigma, 0.0 * unit.kilojoule_per_mole)
    for exception_index in range(nonbonded.getNumExceptions()):
        atom1, atom2, _, sigma, _ = nonbonded.getExceptionParameters(exception_index)
        if atom1 in inactive_indices or atom2 in inactive_indices:
            nonbonded.setExceptionParameters(
                exception_index,
                atom1,
                atom2,
                0.0 * unit.elementary_charge**2,
                sigma,
                0.0 * unit.kilojoule_per_mole,
            )


def build_system(run_dir: Path):
    build = run_dir / "build"
    psf = CharmmPsfFile(str(run_dir / "hybrid.psf"))
    params = CharmmParameterSet(
        str(run_dir.parents[3] / "toppar" / "top_all36_cgenff.rtf"),
        str(run_dir.parents[3] / "toppar" / "par_all36_cgenff.prm"),
        str(build / "core.rtf"),
        str(build / "full_ligand.prm"),
        str(build / "site1_sub1_pres.rtf"),
        str(build / "site1_sub2_pres.rtf"),
    )
    system = psf.createSystem(
        params,
        nonbondedMethod=NoCutoff,
        constraints=HBonds,
        removeCMMotion=True,
    )
    return psf, system


def run_endpoint(
    run_dir: Path,
    endpoint_name: str,
    steps: int,
    report_interval: int,
    seed: int,
) -> dict:
    endpoint = ENDPOINTS[endpoint_name]
    psf, system = build_system(run_dir)
    pdb = PDBFile(str(run_dir / "hybrid.pdb"))
    names = {atom.name: atom.index for atom in pdb.topology.atoms()}
    inactive_indices = {names[name] for name in endpoint["inactive"]}
    deactivate_fragment(system, inactive_indices)

    integrator = LangevinMiddleIntegrator(
        300.0 * unit.kelvin,
        1.0 / unit.picosecond,
        1.0 * unit.femtosecond,
    )
    integrator.setRandomNumberSeed(seed)
    platform = Platform.getPlatformByName("CPU")
    simulation = Simulation(
        psf.topology,
        system,
        integrator,
        platform,
        {"Threads": "1"},
    )
    simulation.context.setPositions(pdb.positions)
    simulation.minimizeEnergy(
        tolerance=10.0 * unit.kilojoule_per_mole / unit.nanometer,
        maxIterations=2000,
    )
    simulation.context.setVelocitiesToTemperature(300.0 * unit.kelvin, seed)

    rows: list[dict[str, float]] = []
    improper_indices = [names[name] for name in endpoint["improper"]]

    def sample(step: int) -> None:
        state = simulation.context.getState(getPositions=True, getEnergy=True)
        positions = np.asarray(
            state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        )
        rows.append(
            {
                "step": step,
                "time_ps": step * 0.001,
                "signed_volume_nm3": signed_volume(positions, names, endpoint),
                "improper_deg": torsion_degrees(positions[improper_indices]),
                "potential_kj_mol": state.getPotentialEnergy().value_in_unit(
                    unit.kilojoule_per_mole
                ),
            }
        )

    sample(0)
    for step in range(report_interval, steps + 1, report_interval):
        simulation.step(report_interval)
        sample(step)

    csv_path = run_dir / f"endpoint_{endpoint_name}.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    volumes = np.asarray([row["signed_volume_nm3"] for row in rows])
    impropers = np.asarray([row["improper_deg"] for row in rows])
    expected_sign = endpoint["expected_volume_sign"]
    sign_flips = int(np.count_nonzero(np.sign(volumes) != expected_sign))
    return {
        "lambda": endpoint["lambda"],
        "steps": steps,
        "duration_ps": steps * 0.001,
        "samples": len(rows),
        "expected_signed_volume_sign": expected_sign,
        "sign_flip_samples": sign_flips,
        "configuration_retained": sign_flips == 0,
        "signed_volume_nm3": {
            "initial": float(volumes[0]),
            "minimum": float(volumes.min()),
            "maximum": float(volumes.max()),
            "mean": float(volumes.mean()),
            "minimum_absolute": float(np.abs(volumes).min()),
        },
        "improper_degrees": {
            "initial": float(impropers[0]),
            "minimum": float(impropers.min()),
            "maximum": float(impropers.max()),
            "mean": float(impropers.mean()),
        },
        "trajectory": csv_path.name,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--steps", type=int, default=1_000_000)
    parser.add_argument("--report-interval", type=int, default=1_000)
    args = parser.parse_args()
    run_dir = args.run_dir.resolve()

    results = {
        "conditions": {
            "temperature_K": 300.0,
            "timestep_fs": 1.0,
            "friction_per_ps": 1.0,
            "nonbonded_method": "NoCutoff (vacuum smoke test)",
            "constraints": "HBonds",
            "platform": "CPU, one thread",
            "improper_force_constant_kcal_mol_rad2": 50.0,
            "improper_target_degrees": 32.60,
        },
        "endpoints": {},
    }
    for index, endpoint_name in enumerate(ENDPOINTS):
        results["endpoints"][endpoint_name] = run_endpoint(
            run_dir,
            endpoint_name,
            args.steps,
            args.report_interval,
            seed=20260919 + index,
        )
    (run_dir / "endpoint_summary.json").write_text(
        json.dumps(results, indent=2) + "\n"
    )
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
