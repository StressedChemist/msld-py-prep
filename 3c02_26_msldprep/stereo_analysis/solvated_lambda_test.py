#!/usr/bin/env python3
"""Solvated two-state lambda/Gibbs test with stereochemistry tracking.

This is a compact validation system: the prepared meso/threitol dual topology
is placed in periodic TIP3P water, alchemical nonbonded interactions are scaled
over an 11-state lambda grid, and a bias-updated Gibbs sampler moves lambda.
Both endpoint-specific C3 signed volumes and improper angles are recorded after
every MD/Gibbs cycle.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from openmm import (
    CustomBondForce,
    CustomNonbondedForce,
    LangevinMiddleIntegrator,
    NonbondedForce,
    Platform,
    Vec3,
    XmlSerializer,
    unit,
)
from openmm.app import (
    CharmmParameterSet,
    CharmmPsfFile,
    CutoffPeriodic,
    DCDReporter,
    ForceField,
    HBonds,
    Modeller,
    PDBFile,
    Simulation,
    Topology,
)


MESO_NAMES = ["O015", "H016", "H017", "C018"]
THREO_NAMES = ["O019", "H020", "H021", "C022"]
R_KJ_MOL_K = 0.00831446261815324


def torsion_degrees(points: np.ndarray) -> float:
    p0, p1, p2, p3 = points
    b0 = -(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return math.degrees(math.atan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))


def signed_volume(
    positions: np.ndarray,
    names: dict[str, int],
    centre_name: str,
    oxygen_name: str,
) -> float:
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


def make_water_box(box_nm: float, ligand_positions_nm: np.ndarray):
    water_forcefield = ForceField("tip3p.xml")
    water = Modeller(Topology(), [])
    water.addSolvent(
        water_forcefield,
        boxSize=Vec3(box_nm, box_nm, box_nm) * unit.nanometer,
        model="tip3p",
        neutralize=False,
    )

    # The empty-box construction has no solute cavity. Remove waters whose
    # oxygen overlaps the centred ligand, then let minimization relax the shell.
    delete_residues = []
    water_positions_nm = np.asarray(
        water.positions.value_in_unit(unit.nanometer)
    )
    for residue in water.topology.residues():
        oxygen = next(atom for atom in residue.atoms() if atom.element.symbol == "O")
        distances = np.linalg.norm(
            ligand_positions_nm - water_positions_nm[oxygen.index], axis=1
        )
        if float(distances.min()) < 0.26:
            delete_residues.append(residue)
    water.delete(delete_residues)
    return water_forcefield, water


def merge_water_system(ligand_system, water_system) -> None:
    offset = ligand_system.getNumParticles()
    for particle_index in range(water_system.getNumParticles()):
        ligand_system.addParticle(water_system.getParticleMass(particle_index))
    for constraint_index in range(water_system.getNumConstraints()):
        atom1, atom2, distance = water_system.getConstraintParameters(constraint_index)
        ligand_system.addConstraint(atom1 + offset, atom2 + offset, distance)

    ligand_nonbonded = next(
        force for force in ligand_system.getForces() if isinstance(force, NonbondedForce)
    )
    water_nonbonded = next(
        force for force in water_system.getForces() if isinstance(force, NonbondedForce)
    )
    for particle_index in range(water_nonbonded.getNumParticles()):
        ligand_nonbonded.addParticle(
            *water_nonbonded.getParticleParameters(particle_index)
        )
    for exception_index in range(water_nonbonded.getNumExceptions()):
        atom1, atom2, chargeprod, sigma, epsilon = (
            water_nonbonded.getExceptionParameters(exception_index)
        )
        ligand_nonbonded.addException(
            atom1 + offset,
            atom2 + offset,
            chargeprod,
            sigma,
            epsilon,
        )


def add_alchemical_nonbonded(system, names: dict[str, int], cutoff_nm: float) -> None:
    original_index, original = next(
        (index, force)
        for index, force in enumerate(system.getForces())
        if isinstance(force, NonbondedForce)
    )
    meso = {names[name] for name in MESO_NAMES}
    threo = {names[name] for name in THREO_NAMES}
    substituent = [0] * system.getNumParticles()
    site = [0] * system.getNumParticles()
    for atom_index in meso:
        substituent[atom_index] = 1
        site[atom_index] = 1
    for atom_index in threo:
        substituent[atom_index] = 2
        site[atom_index] = 1

    expression = (
        "scale*(138.935456*charge1*charge2/rsoft+"
        "4*sqrt(epsilon1*epsilon2)*(sr6*sr6-sr6));"
        "sr6=(sigma/rsoft)^6;"
        "rsoft=(r^6+0.5*(1-scale)*sigma^6)^(1/6);"
        "sigma=0.5*(sigma1+sigma2);"
        "scale=samesite*samesub*w1+(1-samesite)*w1*w2;"
        "samesite=delta(site1-site2);"
        "samesub=delta(sub1-sub2);"
        "w1=delta(sub1)+delta(sub1-1)*lambda_stereo+"
        "delta(sub1-2)*(1-lambda_stereo);"
        "w2=delta(sub2)+delta(sub2-1)*lambda_stereo+"
        "delta(sub2-2)*(1-lambda_stereo)"
    )
    custom = CustomNonbondedForce(expression)
    for parameter in ("charge", "sigma", "epsilon", "sub", "site"):
        custom.addPerParticleParameter(parameter)
    custom.addGlobalParameter("lambda_stereo", 1.0)
    custom.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
    custom.setCutoffDistance(cutoff_nm * unit.nanometer)
    custom.setUseSwitchingFunction(True)
    custom.setSwitchingDistance((cutoff_nm - 0.1) * unit.nanometer)
    custom.setUseLongRangeCorrection(False)

    particles = []
    for atom_index in range(original.getNumParticles()):
        charge, sigma, epsilon = original.getParticleParameters(atom_index)
        particle = [
            charge.value_in_unit(unit.elementary_charge),
            sigma.value_in_unit(unit.nanometer),
            epsilon.value_in_unit(unit.kilojoule_per_mole),
            substituent[atom_index],
            site[atom_index],
        ]
        particles.append(particle)
        custom.addParticle(particle)

    bond_expression = (
        "scale*(138.935456*chargeprod/r+"
        "4*epsilon*((sigma/r)^12-(sigma/r)^6));"
        "scale=samesite*samesub*w1+(1-samesite)*w1*w2;"
        "samesite=delta(site1-site2);"
        "samesub=delta(sub1-sub2);"
        "w1=delta(sub1)+delta(sub1-1)*lambda_stereo+"
        "delta(sub1-2)*(1-lambda_stereo);"
        "w2=delta(sub2)+delta(sub2-1)*lambda_stereo+"
        "delta(sub2-2)*(1-lambda_stereo)"
    )
    one_four = CustomBondForce(bond_expression)
    for parameter in (
        "chargeprod",
        "sigma",
        "epsilon",
        "sub1",
        "sub2",
        "site1",
        "site2",
    ):
        one_four.addPerBondParameter(parameter)
    one_four.addGlobalParameter("lambda_stereo", 1.0)

    exclusions: set[tuple[int, int]] = set()
    for exception_index in range(original.getNumExceptions()):
        atom1, atom2, chargeprod, sigma, epsilon = (
            original.getExceptionParameters(exception_index)
        )
        pair = tuple(sorted((int(atom1), int(atom2))))
        if pair not in exclusions:
            custom.addExclusion(*pair)
            exclusions.add(pair)
        q = chargeprod.value_in_unit(unit.elementary_charge**2)
        e = epsilon.value_in_unit(unit.kilojoule_per_mole)
        cross_alternative = (
            (atom1 in meso and atom2 in threo)
            or (atom1 in threo and atom2 in meso)
        )
        if (q != 0.0 or e != 0.0) and not cross_alternative:
            one_four.addBond(
                atom1,
                atom2,
                [
                    q,
                    sigma.value_in_unit(unit.nanometer),
                    e,
                    substituent[atom1],
                    substituent[atom2],
                    site[atom1],
                    site[atom2],
                ],
            )

    # Opposing alternatives never interact. C018 and C022 start coincident, so
    # explicit exclusions also prevent a zero-times-infinity singularity.
    for atom1 in meso:
        for atom2 in threo:
            pair = tuple(sorted((atom1, atom2)))
            if pair not in exclusions:
                custom.addExclusion(*pair)
                exclusions.add(pair)

    custom.setForceGroup(20)
    one_four.setForceGroup(20)
    system.addForce(custom)
    system.addForce(one_four)
    system.removeForce(original_index)


def build_solvated_system(run_dir: Path, output_dir: Path, box_nm: float):
    build = run_dir / "build"
    psf = CharmmPsfFile(str(run_dir / "hybrid.psf"))
    pdb = PDBFile(str(run_dir / "hybrid.pdb"))
    prep_root = run_dir.parents[3]
    params = CharmmParameterSet(
        str(prep_root / "toppar" / "top_all36_cgenff.rtf"),
        str(prep_root / "toppar" / "par_all36_cgenff.prm"),
        str(build / "core.rtf"),
        str(build / "full_ligand.prm"),
        str(build / "site1_sub1_pres.rtf"),
        str(build / "site1_sub2_pres.rtf"),
    )
    cutoff_nm = min(1.0, box_nm / 2.0 - 0.1)
    psf.setBox(
        box_nm * unit.nanometer,
        box_nm * unit.nanometer,
        box_nm * unit.nanometer,
        90.0,
        90.0,
        90.0,
    )
    system = psf.createSystem(
        params,
        nonbondedMethod=CutoffPeriodic,
        nonbondedCutoff=cutoff_nm * unit.nanometer,
        constraints=HBonds,
        removeCMMotion=True,
    )

    ligand_positions_nm = np.asarray(
        pdb.positions.value_in_unit(unit.nanometer)
    )
    ligand_positions_nm += box_nm / 2.0 - ligand_positions_nm.mean(axis=0)
    water_forcefield, water = make_water_box(box_nm, ligand_positions_nm)
    water_system = water_forcefield.createSystem(
        water.topology,
        nonbondedMethod=CutoffPeriodic,
        nonbondedCutoff=cutoff_nm * unit.nanometer,
        constraints=HBonds,
        rigidWater=True,
        removeCMMotion=False,
    )
    merge_water_system(system, water_system)
    vectors = (
        Vec3(box_nm, 0, 0),
        Vec3(0, box_nm, 0),
        Vec3(0, 0, box_nm),
    ) * unit.nanometer
    system.setDefaultPeriodicBoxVectors(*vectors)

    combined = Modeller(
        psf.topology,
        ligand_positions_nm * unit.nanometer,
    )
    combined.add(water.topology, water.positions)
    combined.topology.setPeriodicBoxVectors(vectors)
    names = {
        atom.name: atom.index
        for atom in combined.topology.atoms()
        if atom.residue.name == "LIG"
    }
    add_alchemical_nonbonded(system, names, cutoff_nm)

    with (output_dir / "solvated_initial.pdb").open("w") as handle:
        PDBFile.writeFile(combined.topology, combined.positions, handle)
    (output_dir / "solvated_system.xml").write_text(
        XmlSerializer.serialize(system)
    )
    return combined, system, names, len(list(water.topology.residues())), cutoff_nm


def endpoint_alternations(states: list[int], last_state: int) -> int:
    extremes = []
    for state in states:
        if state in (0, last_state) and (not extremes or state != extremes[-1]):
            extremes.append(state)
    return max(0, len(extremes) - 1)


def run(args) -> dict:
    run_dir = args.run_dir.resolve()
    output_dir = run_dir / args.output_name
    if output_dir.exists():
        raise RuntimeError(f"Refusing to overwrite existing output: {output_dir}")
    output_dir.mkdir(parents=True)

    modeller, system, names, water_count, cutoff_nm = build_solvated_system(
        run_dir, output_dir, args.box_nm
    )
    print(
        f"built solvated system: {system.getNumParticles()} atoms, "
        f"{water_count} waters, cutoff={cutoff_nm:.2f} nm",
        flush=True,
    )
    integrator = LangevinMiddleIntegrator(
        args.temperature * unit.kelvin,
        1.0 / unit.picosecond,
        args.timestep_fs * unit.femtosecond,
    )
    integrator.setRandomNumberSeed(args.seed)
    platform = Platform.getPlatformByName("CPU")
    simulation = Simulation(
        modeller.topology,
        system,
        integrator,
        platform,
        {"Threads": str(args.threads)},
    )
    simulation.context.setPositions(modeller.positions)
    simulation.context.setParameter("lambda_stereo", 1.0)
    print("minimizing", flush=True)
    simulation.minimizeEnergy(
        tolerance=args.minimization_tolerance
        * unit.kilojoule_per_mole
        / unit.nanometer,
        maxIterations=args.minimization_iterations,
    )
    minimized = simulation.context.getState(
        getPositions=True, getEnergy=True, getForces=True
    )
    minimized_forces = np.asarray(
        minimized.getForces(asNumpy=True).value_in_unit(
            unit.kilojoule_per_mole / unit.nanometer
        )
    )
    force_norms = np.linalg.norm(minimized_forces, axis=1)
    max_force_index = int(force_norms.argmax())
    max_force_atom = list(modeller.topology.atoms())[max_force_index]
    print(
        "minimization complete: "
        f"E={minimized.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole):.3f} "
        f"kJ/mol, max|F|={force_norms[max_force_index]:.3f} kJ/mol/nm "
        f"at {max_force_atom.residue.name}:{max_force_atom.residue.id}:"
        f"{max_force_atom.name}[{max_force_index}]",
        flush=True,
    )
    with (output_dir / "solvated_minimized.pdb").open("w") as handle:
        PDBFile.writeFile(
            modeller.topology,
            minimized.getPositions(),
            handle,
            keepIds=True,
        )
    simulation.context.setVelocitiesToTemperature(
        args.temperature * unit.kelvin, args.seed
    )
    if args.equilibration_steps:
        print(f"equilibrating for {args.equilibration_steps} steps", flush=True)
        simulation.step(args.equilibration_steps)
        print("equilibration complete", flush=True)
    simulation.reporters.append(
        DCDReporter(str(output_dir / "trajectory.dcd"), args.trajectory_interval)
    )

    lambda_values = np.linspace(0.0, 1.0, args.lambda_states)
    state_index = args.lambda_states - 1
    simulation.context.setParameter("lambda_stereo", lambda_values[state_index])
    log_weights = np.zeros(args.lambda_states)
    total_counts = np.zeros(args.lambda_states, dtype=int)
    window_counts = np.zeros(args.lambda_states, dtype=int)
    gamma = args.initial_bias_increment
    beta = 1.0 / (R_KJ_MOL_K * args.temperature)
    rng = np.random.default_rng(args.seed)

    meso_improper = [names[name] for name in ("C018", "C009", "C002", "O015")]
    threo_improper = [names[name] for name in ("C022", "C002", "C009", "O019")]
    rows = []
    state_history = []
    first_energy_span = None
    last_energy_span = None

    for cycle in range(1, args.cycles + 1):
        simulation.step(args.steps_per_cycle)
        energies = np.empty(args.lambda_states)
        for trial_index, trial_lambda in enumerate(lambda_values):
            simulation.context.setParameter("lambda_stereo", float(trial_lambda))
            energies[trial_index] = simulation.context.getState(
                getEnergy=True
            ).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        last_energy_span = float(energies.max() - energies.min())
        if first_energy_span is None:
            first_energy_span = last_energy_span

        log_probability = -beta * energies + log_weights
        log_probability -= log_probability.max()
        probability = np.exp(log_probability)
        probability /= probability.sum()
        state_index = int(rng.choice(args.lambda_states, p=probability))
        selected_lambda = float(lambda_values[state_index])
        simulation.context.setParameter("lambda_stereo", selected_lambda)
        total_counts[state_index] += 1
        window_counts[state_index] += 1
        state_history.append(state_index)

        log_weights[state_index] -= gamma
        log_weights -= log_weights.mean()
        if cycle % args.flatness_interval == 0:
            mean_count = float(window_counts.mean())
            if (
                mean_count > 0
                and float(window_counts.min()) >= args.flatness_fraction * mean_count
                and gamma > args.minimum_bias_increment
            ):
                gamma = max(args.minimum_bias_increment, gamma / 2.0)
                window_counts[:] = 0

        state = simulation.context.getState(getPositions=True)
        positions = np.asarray(
            state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        )
        rows.append(
            {
                "cycle": cycle,
                "time_ps": (
                    args.equilibration_steps + cycle * args.steps_per_cycle
                )
                * args.timestep_fs
                / 1000.0,
                "state": state_index,
                "lambda_meso": selected_lambda,
                "lambda_threo": 1.0 - selected_lambda,
                "selected_probability": float(probability[state_index]),
                "energy_selected_kj_mol": float(energies[state_index]),
                "meso_signed_volume_nm3": signed_volume(
                    positions, names, "C018", "O015"
                ),
                "threo_signed_volume_nm3": signed_volume(
                    positions, names, "C022", "O019"
                ),
                "meso_improper_deg": torsion_degrees(positions[meso_improper]),
                "threo_improper_deg": torsion_degrees(positions[threo_improper]),
                "bias_increment": gamma,
            }
        )
        if cycle % args.progress_interval == 0:
            print(
                f"cycle={cycle} lambda={selected_lambda:.2f} "
                f"visited={np.count_nonzero(total_counts)}/{args.lambda_states} "
                f"transitions={np.count_nonzero(np.diff(state_history))} "
                f"gamma={gamma:.4f}",
                flush=True,
            )

    with (output_dir / "lambda_stereo.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    final_state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
    with (output_dir / "solvated_final.pdb").open("w") as handle:
        PDBFile.writeFile(
            modeller.topology,
            final_state.getPositions(),
            handle,
            keepIds=True,
        )

    meso_volumes = np.asarray([row["meso_signed_volume_nm3"] for row in rows])
    threo_volumes = np.asarray([row["threo_signed_volume_nm3"] for row in rows])
    meso_angles = np.asarray([row["meso_improper_deg"] for row in rows])
    threo_angles = np.asarray([row["threo_improper_deg"] for row in rows])
    transitions = int(np.count_nonzero(np.diff(state_history)))
    summary = {
        "conditions": {
            "water_model": "TIP3P",
            "water_molecules": water_count,
            "box_nm": args.box_nm,
            "cutoff_nm": cutoff_nm,
            "ensemble": "NVT",
            "temperature_K": args.temperature,
            "timestep_fs": args.timestep_fs,
            "equilibration_steps_at_meso_endpoint": args.equilibration_steps,
            "cycles": args.cycles,
            "steps_per_cycle": args.steps_per_cycle,
            "lambda_states": lambda_values.tolist(),
            "sampler": "bias-updated all-state Gibbs",
        },
        "lambda_sampling": {
            "state_counts": total_counts.tolist(),
            "states_visited": int(np.count_nonzero(total_counts)),
            "transitions": transitions,
            "endpoint_alternations": endpoint_alternations(
                state_history, args.lambda_states - 1
            ),
            "final_log_weights": log_weights.tolist(),
            "final_bias_increment": gamma,
            "all_states_visited": bool(np.all(total_counts > 0)),
            "first_cycle_energy_span_kj_mol": first_energy_span,
            "last_cycle_energy_span_kj_mol": last_energy_span,
        },
        "stereochemistry": {
            "meso": {
                "expected_sign": -1,
                "sign_flip_samples": int(np.count_nonzero(meso_volumes >= 0)),
                "configuration_retained": bool(np.all(meso_volumes < 0)),
                "signed_volume_nm3_min": float(meso_volumes.min()),
                "signed_volume_nm3_max": float(meso_volumes.max()),
                "minimum_absolute_volume_nm3": float(np.abs(meso_volumes).min()),
                "improper_deg_min": float(meso_angles.min()),
                "improper_deg_max": float(meso_angles.max()),
            },
            "threo": {
                "expected_sign": 1,
                "sign_flip_samples": int(np.count_nonzero(threo_volumes <= 0)),
                "configuration_retained": bool(np.all(threo_volumes > 0)),
                "signed_volume_nm3_min": float(threo_volumes.min()),
                "signed_volume_nm3_max": float(threo_volumes.max()),
                "minimum_absolute_volume_nm3": float(np.abs(threo_volumes).min()),
                "improper_deg_min": float(threo_angles.min()),
                "improper_deg_max": float(threo_angles.max()),
            },
            "separated_throughout": bool(
                np.all(meso_volumes < 0) and np.all(threo_volumes > 0)
            ),
        },
        "outputs": {
            "timeseries": "lambda_stereo.csv",
            "trajectory": "trajectory.dcd",
            "initial_structure": "solvated_initial.pdb",
            "final_structure": "solvated_final.pdb",
            "system": "solvated_system.xml",
        },
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2), flush=True)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--output-name", default="solvated_lambda")
    parser.add_argument("--box-nm", type=float, default=2.6)
    parser.add_argument("--temperature", type=float, default=300.0)
    parser.add_argument("--timestep-fs", type=float, default=2.0)
    parser.add_argument("--equilibration-steps", type=int, default=5000)
    parser.add_argument("--minimization-iterations", type=int, default=250)
    parser.add_argument("--minimization-tolerance", type=float, default=100.0)
    parser.add_argument("--cycles", type=int, default=5000)
    parser.add_argument("--steps-per-cycle", type=int, default=100)
    parser.add_argument("--lambda-states", type=int, default=11)
    parser.add_argument("--initial-bias-increment", type=float, default=0.5)
    parser.add_argument("--minimum-bias-increment", type=float, default=0.01)
    parser.add_argument("--flatness-interval", type=int, default=250)
    parser.add_argument("--flatness-fraction", type=float, default=0.6)
    parser.add_argument("--trajectory-interval", type=int, default=500)
    parser.add_argument("--progress-interval", type=int, default=100)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--seed", type=int, default=20260919)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
