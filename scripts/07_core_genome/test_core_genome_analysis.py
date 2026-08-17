#!/usr/bin/env python3
"""Small static unit tests for accession and QC helpers."""

import importlib.util
import unittest
from pathlib import Path

import pandas as pd
import numpy as np

MODULE_PATH = Path(__file__).with_name("core_genome_analysis.py")
SPEC = importlib.util.spec_from_file_location("core_genome_analysis", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class CoreHelpersTest(unittest.TestCase):
    def test_versioned_accession_is_preserved(self):
        self.assertEqual(MODULE.extract_accession("reproductive__GCA_029011255.1.gff"), "GCA_029011255.1")

    def test_unversioned_accession_is_rejected(self):
        with self.assertRaises(ValueError):
            MODULE.extract_accession("GCA_029011255.gff")

    def test_iqr_flags_do_not_remove_rows(self):
        frame = pd.DataFrame({"metric": [1, 1, 1, 100]})
        flagged = MODULE.add_iqr_flags(frame, ["metric"])
        self.assertEqual(len(flagged), 4)
        self.assertIn("outside_1.5_IQR", flagged.iloc[-1].Outlier_flags)

    def test_all_tied_minimum_neighbours_are_retained(self):
        distances = pd.Series([0, 3, 3, 8], index=["A", "B", "C", "D"])
        minimum, tied = MODULE.tied_minimum_neighbours(distances, "A")
        self.assertEqual(minimum, 3)
        self.assertEqual(tied, ["B", "C"])

    def test_distance_matrix_rejects_asymmetry(self):
        matrix = pd.DataFrame([[0, 1], [2, 0]], index=["A", "B"], columns=["A", "B"])
        with self.assertRaisesRegex(ValueError, "symmetric"):
            MODULE.validate_distance_matrix(matrix, {"A", "B"})

    def test_distance_matrix_rejects_missing_genome(self):
        matrix = pd.DataFrame(np.zeros((2, 2)), index=["A", "B"], columns=["A", "B"])
        with self.assertRaisesRegex(ValueError, "canonical accession"):
            MODULE.validate_distance_matrix(matrix, {"A", "B", "C"})


if __name__ == "__main__":
    unittest.main()
