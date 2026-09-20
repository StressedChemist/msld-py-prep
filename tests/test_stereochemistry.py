#!/usr/bin/env python3
"""Regression tests for chirality-aware MCS preparation."""

from __future__ import annotations

import unittest
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdFMCS

from msld_mcs_rdecomp import assign_and_validate_stereochemistry


class StereochemistryTests(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
