"""Synthetic tests for shared-shift projection and own-target mapping.

Expected values are literals in this file. Tests do not read the Norman H5AD
and do not refit models.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from modeling_diagnostics import (
    map_own_target_indices,
    own_target_value,
    project_residual_onto_shift,
    select_example_pairs,
    top_cosine_pairs,
)


class TestProjectionOntoSharedShift(unittest.TestCase):
    def test_parallel_residual_has_unit_squared_norm_fraction(self):
        direction = np.array([1.0, 0.0, 0.0])
        residual = np.array([3.0, 0.0, 0.0])
        out = project_residual_onto_shift(residual, direction)
        self.assertAlmostEqual(out["coefficient"], 3.0)
        self.assertAlmostEqual(out["squared_norm_fraction"], 1.0)
        self.assertAlmostEqual(out["cosine_to_neg_direction"], -1.0)
        self.assertAlmostEqual(out["orthogonal_rms"], 0.0)

    def test_orthogonal_residual_has_zero_squared_norm_fraction(self):
        direction = np.array([1.0, 0.0])
        residual = np.array([0.0, 2.0])
        out = project_residual_onto_shift(residual, direction)
        self.assertAlmostEqual(out["coefficient"], 0.0)
        self.assertAlmostEqual(out["squared_norm_fraction"], 0.0)
        self.assertAlmostEqual(out["cosine_to_neg_direction"], 0.0)
        self.assertAlmostEqual(out["orthogonal_rms"], 2.0 / np.sqrt(2.0))

    def test_zero_residual_or_direction_is_nan(self):
        zero = np.array([0.0, 0.0, 0.0])
        other = np.array([1.0, 2.0, 3.0])
        for residual, direction in ((zero, other), (other, zero), (zero, zero)):
            out = project_residual_onto_shift(residual, direction)
            self.assertTrue(np.isnan(out["coefficient"]))
            self.assertTrue(np.isnan(out["squared_norm_fraction"]))
            self.assertTrue(np.isnan(out["cosine_to_neg_direction"]))
            self.assertTrue(np.isnan(out["orthogonal_rms"]))


class TestOwnTargetMapping(unittest.TestCase):
    def test_exact_gene_name_maps_to_feature_index(self):
        mapping, missing = map_own_target_indices(
            ("CEBPA", "HBZ", "MISSING"),
            ("LST1", "CEBPA", "HBZ"),
        )
        self.assertEqual(mapping["CEBPA"], 1)
        self.assertEqual(mapping["HBZ"], 2)
        self.assertIsNone(mapping["MISSING"])
        self.assertEqual(missing, ["MISSING"])

    def test_own_target_reads_mapped_coordinate_not_first_gene(self):
        effects = np.array([0.1, 0.8, -0.2])
        mapping, missing = map_own_target_indices(("B",), ("A", "B", "C"))
        self.assertEqual(missing, [])
        self.assertAlmostEqual(own_target_value(effects, mapping["B"]), 0.8)
        self.assertNotAlmostEqual(own_target_value(effects, mapping["B"]), 0.1)

    def test_duplicate_gene_name_is_rejected(self):
        mapping, missing = map_own_target_indices(("A",), ("A", "A"))
        self.assertIsNone(mapping["A"])
        self.assertEqual(missing, ["A"])


class TestNeighborAndExampleSelection(unittest.TestCase):
    def test_top_cosine_pairs_exclude_diagonal(self):
        signatures = np.array(
            [
                [1.0, 0.0, 0.0],
                [0.99, 0.1, 0.0],
                [0.0, 0.0, 1.0],
            ]
        )
        pairs = top_cosine_pairs(signatures, ("A", "B", "C"), n_pairs=1)
        self.assertEqual(len(pairs), 1)
        self.assertEqual(set(pairs[0]["pair"]), {"A", "B"})
        self.assertLess(pairs[0]["cosine"], 1.0)
        self.assertGreater(pairs[0]["cosine"], 0.9)

    def test_example_pairs_prefer_distinct_roles(self):
        chosen = select_example_pairs(
            pair_keys=["P1", "P2", "P3", "P4", "P5"],
            residual_rms=np.array([0.40, 0.20, 0.30, 0.10, 0.35]),
            observed_rms=np.array([0.50, 0.50, 0.05, 0.60, 0.55]),
            relative_residual=np.array([0.80, 0.40, 6.00, 0.16, 0.64]),
            top3_min_sign=np.array([0.50, 0.50, 0.50, 1.00, 1.00]),
            min_sign=0.875,
        )
        self.assertEqual([row["role"] for row in chosen], [
            "largest_canonical_residual",
            "lowest_relative_residual_above_median_observed_rms",
            "largest_lane_sign_supported_residual",
        ])
        self.assertEqual(chosen[0]["canonical_perturbation"], "P1")
        self.assertEqual(chosen[1]["canonical_perturbation"], "P4")
        self.assertEqual(chosen[2]["canonical_perturbation"], "P5")
        self.assertEqual(len({row["canonical_perturbation"] for row in chosen}), 3)


if __name__ == "__main__":
    unittest.main()
