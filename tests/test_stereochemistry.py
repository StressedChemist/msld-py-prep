#!/usr/bin/env python3
"""Regression tests for chirality-aware MCS preparation."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdFMCS

from msld_mcs_rdecomp import (
    assign_and_validate_stereochemistry,
    write_stereochemistry_audit,
)


class StereochemistryTests(unittest.TestCase):
    @staticmethod
    def _erythritol_pair():
        inputs = (
            Path(__file__).resolve().parents[1]
            / "3c02_26_msldprep"
            / "stereo_analysis"
            / "inputs"
        )
        stems = ("meso_erythritol", "threitol")
        molecules = [
            Chem.SDMolSupplier(str(inputs / f"{stem}.sdf"), removeHs=False)[0]
            for stem in stems
        ]
        centres = [
            assign_and_validate_stereochemistry(molecule, stem)
            for molecule, stem in zip(molecules, stems)
        ]
        atom_names = []
        atom_types = []
        for stem in stems:
            names = []
            types = []
            for line in (inputs / f"{stem}.rtf").read_text().splitlines():
                fields = line.split("!", 1)[0].split()
                if fields and fields[0] == "ATOM" and not fields[1].startswith("LP"):
                    names.append(fields[1])
                    types.append(fields[2])
            atom_names.append(names)
            atom_types.append(types)
        mcs = rdFMCS.FindMCS(molecules, matchChiralTag=True)
        mcs_indices = [
            molecule.GetSubstructMatch(mcs.queryMol, useChirality=True)
            for molecule in molecules
        ]
        topology_names = [str(inputs / f"{stem}.rtf") for stem in stems]
        parameter_names = [str(inputs / f"{stem}.prm") for stem in stems]
        return (
            stems,
            molecules,
            centres,
            atom_names,
            atom_types,
            mcs_indices,
            topology_names,
            parameter_names,
        )

    def test_shipped_entry_points_use_chirality_aware_mcs(self) -> None:
        repository = Path(__file__).resolve().parents[1]
        entry_points = (
            repository / "msld_py_prep.py",
            repository / "CRN_Plugin" / "msld_py_prep" / "msld_py_prep.py",
            repository / "CRN_Plugin" / "__init__.py",
        )

        for entry_point in entry_points:
            source = entry_point.read_text()
            with self.subTest(entry_point=entry_point):
                self.assertIn("MCSS_RDecomp", source)
                self.assertNotIn("msld_mcs.MsldMCS", source)

        canonical_mcs = (repository / "msld_mcs_rdecomp.py").read_bytes()
        copied_mcs_modules = (
            repository / "3c02_26_msldprep" / "msld_mcs_rdecomp.py",
            repository / "CRN_Plugin" / "msld_py_prep" / "msld_mcs_rdecomp.py",
        )
        for module in copied_mcs_modules:
            with self.subTest(module=module):
                self.assertEqual(module.read_bytes(), canonical_mcs)

    def test_unassigned_tetrahedral_centre_fails_closed(self) -> None:
        molecule = Chem.MolFromSmiles("CCC(O)C")

        with self.assertRaisesRegex(ValueError, "unassigned stereocentres"):
            assign_and_validate_stereochemistry(molecule, "unassigned test")

    def test_diastereomers_do_not_collapse_to_a_full_mcs(self) -> None:
        first = Chem.MolFromSmiles("C[C@H](Cl)[C@H](Cl)C")
        second = Chem.MolFromSmiles("C[C@H](Cl)[C@@H](Cl)C")
        for index, molecule in enumerate((first, second), start=1):
            assign_and_validate_stereochemistry(molecule, f"stereoisomer {index}")

        nonchiral = rdFMCS.FindMCS([first, second], matchChiralTag=False)
        chiral = rdFMCS.FindMCS([first, second], matchChiralTag=True)

        self.assertEqual(nonchiral.numAtoms, first.GetNumAtoms())
        self.assertLess(chiral.numAtoms, first.GetNumAtoms())
        self.assertTrue(first.GetSubstructMatch(chiral.queryMol, useChirality=True))
        self.assertTrue(second.GetSubstructMatch(chiral.queryMol, useChirality=True))

    def test_audit_classifies_pair_and_finds_endpoint_impropers(self) -> None:
        (
            stems,
            molecules,
            centres,
            names,
            types,
            indices,
            topologies,
            parameters,
        ) = self._erythritol_pair()
        with tempfile.TemporaryDirectory() as directory:
            audit, pairs = write_stereochemistry_audit(
                molecules,
                stems,
                names,
                centres,
                indices,
                topologies,
                atom_types=types,
                parameter_names=parameters,
                audit_output=str(Path(directory) / "audit.csv"),
                pairs_output=str(Path(directory) / "pairs.csv"),
            )

        self.assertEqual(pairs[0]["relationship"], "diastereomers")
        self.assertEqual(pairs[0]["retained_centres"], "C2")
        self.assertEqual(pairs[0]["inverted_centres_molecule_1"], "C3")
        changed = [row for row in audit if row["configuration_changes"]]
        self.assertEqual(len(changed), 2)
        self.assertTrue(all(row["status"] == "ready" for row in changed))

    def test_audit_fails_when_changing_centre_has_no_improper(self) -> None:
        (
            stems,
            molecules,
            centres,
            names,
            types,
            indices,
            _,
            parameters,
        ) = self._erythritol_pair()
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            topologies = []
            for stem in stems:
                topology = directory / f"{stem}.rtf"
                topology.write_text("* deliberately lacks an improper\n")
                topologies.append(str(topology))

            with self.assertRaisesRegex(ValueError, "no center-first IMPR record"):
                write_stereochemistry_audit(
                    molecules,
                    stems,
                    names,
                    centres,
                    indices,
                    topologies,
                    atom_types=types,
                    parameter_names=parameters,
                    audit_output=str(directory / "audit.csv"),
                    pairs_output=str(directory / "pairs.csv"),
                )

            self.assertIn("missing_centered_improper", (directory / "audit.csv").read_text())

    def test_audit_labels_complete_inversion_as_enantiomer_candidate(self) -> None:
        molecules = [
            Chem.MolFromSmiles("C[C@H](O)C(=O)O"),
            Chem.MolFromSmiles("C[C@@H](O)C(=O)O"),
        ]
        stems = ("lactate_s", "lactate_r")
        centres = [
            assign_and_validate_stereochemistry(molecule, stem)
            for molecule, stem in zip(molecules, stems)
        ]
        names = [
            [f"A{index}" for index in range(molecule.GetNumAtoms())]
            for molecule in molecules
        ]
        types = [
            [f"T{index}" for index in range(molecule.GetNumAtoms())]
            for molecule in molecules
        ]
        mcs = rdFMCS.FindMCS(molecules, matchChiralTag=True)
        indices = [
            molecule.GetSubstructMatch(mcs.queryMol, useChirality=True)
            for molecule in molecules
        ]

        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            topologies = []
            parameters = []
            for stem in stems:
                topology = directory / f"{stem}.rtf"
                topology.write_text("IMPR A1 A0 A2 A3\n")
                topologies.append(str(topology))
                parameter = directory / f"{stem}.prm"
                parameter.write_text(
                    "IMPROPERS\nT1 T0 T2 T3 50.0 0 35.0\nEND\n"
                )
                parameters.append(str(parameter))
            _, pairs = write_stereochemistry_audit(
                molecules,
                stems,
                names,
                centres,
                indices,
                topologies,
                atom_types=types,
                parameter_names=parameters,
                audit_output=str(directory / "audit.csv"),
                pairs_output=str(directory / "pairs.csv"),
            )

        self.assertEqual(pairs[0]["relationship"], "enantiomer_candidate")

    def test_audit_rejects_planar_improper_target(self) -> None:
        (
            stems,
            molecules,
            centres,
            names,
            types,
            indices,
            _,
            _,
        ) = self._erythritol_pair()
        improper_names = (
            ("C3", "C4", "C2", "O3"),
            ("C3", "C2", "C4", "O3"),
        )
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            topologies = []
            parameters = []
            for molecule_index, stem in enumerate(stems):
                improper = improper_names[molecule_index]
                topology = directory / f"{stem}.rtf"
                topology.write_text("IMPR " + " ".join(improper) + "\n")
                topologies.append(str(topology))
                type_by_name = dict(zip(names[molecule_index], types[molecule_index]))
                improper_types = [type_by_name[name] for name in improper]
                parameter = directory / f"{stem}.prm"
                parameter.write_text(
                    "IMPROPERS\n"
                    + " ".join(improper_types)
                    + " 50.0 0 0.0\nEND\n"
                )
                parameters.append(str(parameter))

            with self.assertRaisesRegex(ValueError, "non-planar target"):
                write_stereochemistry_audit(
                    molecules,
                    stems,
                    names,
                    centres,
                    indices,
                    topologies,
                    atom_types=types,
                    parameter_names=parameters,
                    audit_output=str(directory / "audit.csv"),
                    pairs_output=str(directory / "pairs.csv"),
                )


if __name__ == "__main__":
    unittest.main()
