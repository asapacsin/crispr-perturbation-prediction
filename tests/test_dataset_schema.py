"""Synthetic tests for condition canonicalization, split leakage, and HDF5 categoricals.

Expected values are literals in this file. They are not copied from the
implementation under test.
"""

from __future__ import annotations

import math
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from dataset_schema_lib import (
    InvalidConditionLabel,
    canonicalize_condition,
    combination_split_errors,
    gene_partition_errors,
    parse_condition_tokens,
    target_gene_split_errors,
    target_set,
)
from validate_dataset_schema import SchemaValidationError, _read_anndata_vector


class TestOrientationEquivalence(unittest.TestCase):
    def test_single_target_ctrl_suffix_and_prefix_match(self):
        self.assertEqual(canonicalize_condition("TSC22D1+ctrl"), "TSC22D1")
        self.assertEqual(canonicalize_condition("ctrl+TSC22D1"), "TSC22D1")
        self.assertEqual(
            canonicalize_condition("TSC22D1+ctrl"),
            canonicalize_condition("ctrl+TSC22D1"),
        )

    def test_double_target_order_is_sorted_and_equivalent(self):
        self.assertEqual(canonicalize_condition("KLF1+MAP2K6"), "KLF1+MAP2K6")
        self.assertEqual(canonicalize_condition("MAP2K6+KLF1"), "KLF1+MAP2K6")
        self.assertEqual(
            target_set("KLF1+MAP2K6"),
            frozenset({"KLF1", "MAP2K6"}),
        )
        self.assertEqual(target_set("KLF1+MAP2K6"), target_set("MAP2K6+KLF1"))

    def test_source_RHOXF2_alias_maps_only_when_requested(self):
        self.assertEqual(canonicalize_condition("RHOXF2+ctrl"), "RHOXF2")
        self.assertEqual(
            canonicalize_condition("RHOXF2+ctrl", apply_source_alias=True),
            "RHOXF2BB",
        )
        self.assertEqual(
            canonicalize_condition("RHOXF2+SET", apply_source_alias=True),
            "RHOXF2BB+SET",
        )


class TestCtrlHandling(unittest.TestCase):
    def test_bare_ctrl_is_control(self):
        self.assertEqual(canonicalize_condition("ctrl"), "ctrl")
        self.assertEqual(parse_condition_tokens("ctrl"), [])
        self.assertEqual(target_set("ctrl"), frozenset())

    def test_ctrl_plus_ctrl_is_still_control(self):
        self.assertEqual(canonicalize_condition("ctrl+ctrl"), "ctrl")

    def test_ctrl_token_is_dropped_not_counted_as_target(self):
        self.assertEqual(parse_condition_tokens("CEBPE+ctrl"), ["CEBPE"])
        self.assertNotIn("ctrl", target_set("ctrl+CEBPE"))

    def test_construct_halves_are_not_four_perturbations(self):
        self.assertEqual(len(target_set("KLF1+MAP2K6")), 2)
        self.assertEqual(canonicalize_condition("NegCtrl10+NegCtrl0"), "NegCtrl0+NegCtrl10")


class TestNullAndInvalidLabels(unittest.TestCase):
    def test_none_is_invalid(self):
        with self.assertRaises(InvalidConditionLabel):
            canonicalize_condition(None)

    def test_nan_is_invalid(self):
        with self.assertRaises(InvalidConditionLabel):
            canonicalize_condition(float("nan"))
        self.assertTrue(math.isnan(float("nan")))

    def test_pandas_na_is_invalid_without_ambiguous_boolean(self):
        with self.assertRaises(InvalidConditionLabel):
            canonicalize_condition(pd.NA)

    def test_empty_and_whitespace_are_invalid(self):
        for label in ("", "   ", "\t"):
            with self.assertRaises(InvalidConditionLabel):
                canonicalize_condition(label)

    def test_empty_tokens_are_invalid(self):
        for label in ("GENE+", "+GENE", "++", "KLF1++MAP2K6"):
            with self.assertRaises(InvalidConditionLabel):
                canonicalize_condition(label)

    def test_non_string_is_invalid(self):
        with self.assertRaises(InvalidConditionLabel):
            canonicalize_condition(123)


class TestCombinationSplitLeakage(unittest.TestCase):
    def test_same_unordered_pair_in_train_and_test_is_leakage(self):
        rows = (
            ("KLF1+MAP2K6", "train"),
            ("MAP2K6+KLF1", "test"),
        )
        errors = combination_split_errors(rows)
        self.assertTrue(errors)
        self.assertTrue(any("KLF1+MAP2K6" in err for err in errors))

    def test_same_single_target_set_in_two_splits_is_leakage(self):
        rows = (
            ("TSC22D1+ctrl", "train"),
            ("ctrl+TSC22D1", "validation"),
        )
        errors = combination_split_errors(rows)
        self.assertTrue(errors)

    def test_disjoint_pairs_and_shared_controls_are_allowed(self):
        rows = (
            ("KLF1+MAP2K6", "train"),
            ("FOXA1+FOXL2", "test"),
            ("CEBPA+CEBPB", "validation"),
            ("ctrl", "train"),
            ("ctrl", "validation"),
            ("ctrl", "test"),
        )
        self.assertEqual(combination_split_errors(rows), [])

    def test_single_and_pair_sharing_a_gene_are_not_combination_leaks(self):
        rows = (
            ("KLF1+ctrl", "train"),
            ("MAP2K6+ctrl", "train"),
            ("KLF1+MAP2K6", "test"),
        )
        self.assertEqual(combination_split_errors(rows), [])

    def test_unknown_split_is_rejected_even_for_ctrl(self):
        errors = combination_split_errors((("ctrl", "holdout"),))
        self.assertTrue(errors)
        self.assertTrue(any("unknown split" in err and "holdout" in err for err in errors))


class TestTargetGeneSplitLeakage(unittest.TestCase):
    def test_train_cannot_contain_held_out_test_gene(self):
        rows = (("FOXA1+KLF1", "train"),)
        gene_split = {"FOXA1": "test", "KLF1": "train"}
        errors = target_gene_split_errors(rows, gene_split)
        self.assertTrue(errors)
        self.assertTrue(any("FOXA1" in err for err in errors))

    def test_validation_cannot_contain_test_gene(self):
        rows = (("FOXA1+ctrl", "validation"),)
        gene_split = {"FOXA1": "test"}
        errors = target_gene_split_errors(rows, gene_split)
        self.assertTrue(errors)

    def test_test_cannot_contain_validation_gene(self):
        rows = (("FOXF1+ctrl", "test"),)
        gene_split = {"FOXF1": "validation"}
        errors = target_gene_split_errors(rows, gene_split)
        self.assertTrue(errors)

    def test_val_test_bridge_assigned_to_test_is_leakage(self):
        rows = (("FOXA1+FOXF1", "test"),)
        gene_split = {"FOXA1": "test", "FOXF1": "validation"}
        errors = target_gene_split_errors(rows, gene_split)
        self.assertTrue(errors)
        self.assertTrue(any("bridge" in err.lower() for err in errors))

    def test_explicit_bridge_exclusion_is_allowed(self):
        rows = (("FOXA1+FOXF1", "excluded_val_test_bridge"),)
        gene_split = {"FOXA1": "test", "FOXF1": "validation"}
        self.assertEqual(target_gene_split_errors(rows, gene_split), [])

    def test_control_cannot_be_excluded_as_gene_bridge(self):
        errors = target_gene_split_errors([("ctrl", "excluded_val_test_bridge")], {"A": "train"})
        self.assertTrue(errors)

    def test_bridge_must_contain_validation_and_test_genes(self):
        rows = (("FOXA1+KLF1", "excluded_val_test_bridge"),)
        gene_split = {"FOXA1": "test", "KLF1": "train"}
        errors = target_gene_split_errors(rows, gene_split)
        self.assertTrue(errors)
        self.assertTrue(any("must contain validation AND test" in err for err in errors))

    def test_train_gene_may_appear_as_partner_in_test_combination(self):
        rows = (("FOXA1+KLF1", "test"),)
        gene_split = {"FOXA1": "test", "KLF1": "train"}
        self.assertEqual(target_gene_split_errors(rows, gene_split), [])

    def test_val_or_test_noncontrol_must_contain_heldout_gene(self):
        gene_split = {"FOXA1": "test", "KLF1": "train", "MAP2K6": "train"}
        test_errors = target_gene_split_errors((("KLF1+MAP2K6", "test"),), gene_split)
        val_errors = target_gene_split_errors((("KLF1+ctrl", "validation"),), gene_split)
        self.assertTrue(any("no test-held-out gene" in err for err in test_errors))
        self.assertTrue(any("no validation-held-out gene" in err for err in val_errors))

    def test_unknown_target_gene_is_rejected_in_gene_partition(self):
        rows = (("FOOBAR+KLF1", "train"),)
        gene_split = {"KLF1": "train"}
        errors = target_gene_split_errors(rows, gene_split)
        self.assertTrue(any("FOOBAR" in err and "not in gene partition" in err for err in errors))

    def test_unknown_split_is_rejected_even_for_ctrl(self):
        errors = target_gene_split_errors((("ctrl", "holdout"),), {"FOXA1": "test"})
        self.assertTrue(any("unknown gene-split assignment" in err for err in errors))

    def test_controls_are_ignored_for_gene_leakage(self):
        rows = (("ctrl", "train"), ("ctrl", "test"))
        gene_split = {"FOXA1": "test"}
        self.assertEqual(target_gene_split_errors(rows, gene_split), [])


class TestGenePartition(unittest.TestCase):
    def test_disjoint_labeled_partition_is_ok(self):
        self.assertEqual(
            gene_partition_errors({"KLF1": "train", "FOXA1": "test", "FOXF1": "validation"}),
            [],
        )

    def test_unknown_partition_label_is_rejected(self):
        errors = gene_partition_errors({"KLF1": "holdout"})
        self.assertTrue(any("unknown gene-partition label" in err for err in errors))

    def test_overlapping_assignment_is_rejected(self):
        class TwoMaps(dict):
            def items(self):
                return (("KLF1", "train"), ("KLF1", "test"))

        errors = gene_partition_errors(TwoMaps())
        self.assertTrue(any("overlap" in err for err in errors))


class TestHdf5CategoricalReader(unittest.TestCase):
    def test_legacy_categories_dataset(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "legacy.h5"
            with h5py.File(path, "w") as handle:
                obs = handle.create_group("obs")
                cats = obs.create_group("__categories")
                cats.create_dataset(
                    "condition",
                    data=np.array(["ctrl", "KLF1+ctrl"], dtype="S20"),
                )
                obs.create_dataset("condition", data=np.array([1, 0, 1], dtype="int8"))
            with h5py.File(path, "r") as handle:
                values = _read_anndata_vector(handle["obs"], "condition")
        self.assertEqual(values, ["KLF1+ctrl", "ctrl", "KLF1+ctrl"])

    def test_modern_categorical_group(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "modern.h5"
            with h5py.File(path, "w") as handle:
                obs = handle.create_group("obs")
                condition = obs.create_group("condition")
                condition.create_dataset(
                    "categories",
                    data=np.array(["ctrl", "KLF1+ctrl"], dtype="S20"),
                )
                condition.create_dataset("codes", data=np.array([0, 1], dtype="int8"))
            with h5py.File(path, "r") as handle:
                values = _read_anndata_vector(handle["obs"], "condition")
        self.assertEqual(values, ["ctrl", "KLF1+ctrl"])

    def test_out_of_bounds_code_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "oob.h5"
            with h5py.File(path, "w") as handle:
                obs = handle.create_group("obs")
                cats = obs.create_group("__categories")
                cats.create_dataset("condition", data=np.array(["ctrl"], dtype="S8"))
                obs.create_dataset("condition", data=np.array([2], dtype="int8"))
            with h5py.File(path, "r") as handle:
                with self.assertRaises(SchemaValidationError) as raised:
                    _read_anndata_vector(handle["obs"], "condition")
        self.assertIn("out of bounds", str(raised.exception))

    def test_null_code_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "null.h5"
            with h5py.File(path, "w") as handle:
                obs = handle.create_group("obs")
                cats = obs.create_group("__categories")
                cats.create_dataset("condition", data=np.array(["ctrl"], dtype="S8"))
                obs.create_dataset("condition", data=np.array([-1], dtype="int8"))
            with h5py.File(path, "r") as handle:
                with self.assertRaises(SchemaValidationError) as raised:
                    _read_anndata_vector(handle["obs"], "condition")
        self.assertIn("null categorical code", str(raised.exception))


if __name__ == "__main__":
    unittest.main()
