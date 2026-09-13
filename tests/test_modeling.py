"""Synthetic tests for lane-balanced effects, interaction residuals, and ridge protocol.

Expected values are literals computed independently in this file. They are not
copied from the implementation under test. Tests do not read the Norman H5AD.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from dataset_schema_lib import canonicalize_condition
from modeling_lib import (
    SIGNIFICANCE_NOT_ESTIMABLE,
    additive_delta,
    cosine_similarity,
    delta_pearson,
    fit_ridge_centered,
    group_orientation_labels,
    lane_balanced_mean,
    lane_mean_sampling_variance,
    leave_one_lane_out_rms,
    matching_single_delta,
    multi_hot_rows,
    pair_interaction_lane,
    predict_ridge,
    run_ridge_selection,
    se_of_averaged_lane_contrasts,
    single_effect_lane,
    skill_score,
    summarize_lane_contrasts,
    vector_rms,
    weighted_cell_mean,
    zero_safe_ratio,
)


class TestAdditiveKnownVectorsGiveZeroInteraction(unittest.TestCase):
    def test_additive_pair_residual_is_exactly_zero(self):
        ctrl = np.array([1.0, 2.0, 3.0])
        gene_a = np.array([2.0, 3.0, 4.0])
        gene_b = np.array([3.0, 4.0, 5.0])
        pair = np.array([4.0, 5.0, 6.0])
        residual = pair_interaction_lane(pair, gene_a, gene_b, ctrl)
        np.testing.assert_array_equal(residual, np.array([0.0, 0.0, 0.0]))
        np.testing.assert_array_equal(
            additive_delta(gene_a - ctrl, gene_b - ctrl),
            pair - ctrl,
        )


class TestKnownInteractionResidualRecovered(unittest.TestCase):
    def test_one_gene_excess_is_recovered(self):
        ctrl = np.array([1.0, 2.0, 3.0])
        gene_a = np.array([2.0, 3.0, 4.0])
        gene_b = np.array([3.0, 4.0, 5.0])
        pair = np.array([5.0, 5.0, 6.0])
        residual = pair_interaction_lane(pair, gene_a, gene_b, ctrl)
        np.testing.assert_array_equal(residual, np.array([1.0, 0.0, 0.0]))


class TestLaneBalanceVersusUnequalCellWeight(unittest.TestCase):
    def test_one_large_lane_and_nine_zero_cells_disagree(self):
        lane_means = np.array([[10.0], [0.0]], dtype=np.float64)
        cell_values = np.array(
            [10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            dtype=np.float64,
        )
        balanced = lane_balanced_mean(lane_means)
        weighted = weighted_cell_mean(cell_values.reshape(-1, 1))
        np.testing.assert_allclose(balanced, np.array([5.0]))
        np.testing.assert_allclose(weighted, np.array([1.0]))
        self.assertGreater(float(np.abs(balanced - weighted).item()), 3.9)

    def test_single_effect_uses_lane_means_not_cell_pool(self):
        single_lane = np.array([[4.0], [0.0]])
        ctrl_lane = np.array([[1.0], [1.0]])
        effect = single_effect_lane(single_lane, ctrl_lane)
        np.testing.assert_allclose(lane_balanced_mean(effect), np.array([1.0]))
        pooled_single = (4.0 + 0.0 * 9.0) / 10.0
        pooled_ctrl = 1.0
        self.assertNotAlmostEqual(pooled_single - pooled_ctrl, 1.0)


class TestVarianceContrastIncludesControlOnce(unittest.TestCase):
    def test_two_lane_hand_calculated_standard_error(self):
        # Lane 1 single {4,6}: mean 5, sample var 2, var(mean)=1
        # Lane 1 ctrl   {1,3}: mean 2, sample var 2, var(mean)=1
        # Lane 2 single {7,9}: mean 8, sample var 2, var(mean)=1
        # Lane 2 ctrl   {2,4}: mean 3, sample var 2, var(mean)=1
        # Contrast vars = 2 and 2. Average contrast = 4.
        # SE = sqrt((2+2) / 4) = 1
        lane_contrast_vars = np.array([[2.0], [2.0]])
        se = se_of_averaged_lane_contrasts(lane_contrast_vars)
        np.testing.assert_allclose(se, np.array([1.0]))

    def test_eight_identical_contrast_variances_use_divisor_64(self):
        # 8 lanes, each contrast var 4. SE = sqrt(32/64) = sqrt(1/2)
        lane_contrast_vars = np.full((8, 1), 4.0)
        se = se_of_averaged_lane_contrasts(lane_contrast_vars)
        np.testing.assert_allclose(se, np.array([np.sqrt(0.5)]))

    def test_pair_contrast_adds_control_variance_once(self):
        var_pair = np.array([[1.0], [1.0]])
        var_a = np.array([[1.0], [1.0]])
        var_b = np.array([[1.0], [1.0]])
        var_ctrl = np.array([[1.0], [1.0]])
        summed_once = var_pair + var_a + var_b + var_ctrl
        np.testing.assert_allclose(summed_once, np.array([[4.0], [4.0]]))
        se_once = se_of_averaged_lane_contrasts(summed_once)
        se_if_control_twice = se_of_averaged_lane_contrasts(summed_once + var_ctrl)
        np.testing.assert_allclose(se_once, np.array([np.sqrt(8.0 / 4.0)]))
        self.assertGreater(float(se_if_control_twice[0]), float(se_once[0]))

    def test_sampling_variance_is_sample_var_over_n(self):
        sample_var = np.array([[2.0, 8.0], [3.0, 12.0]])
        n = np.array([2, 3])
        got = lane_mean_sampling_variance(sample_var, n)
        np.testing.assert_allclose(got, np.array([[1.0, 4.0], [1.0, 4.0]]))

    def test_n_less_than_two_sets_se_na_and_flag(self):
        contrasts = np.ones((8, 2))
        ns_ok = [np.full(8, 5), np.full(8, 5)]
        ns_bad = [np.array([5, 5, 5, 1, 5, 5, 5, 5]), np.full(8, 5)]
        vars_ = [np.ones((8, 2)), np.ones((8, 2))]
        ok = summarize_lane_contrasts(contrasts, ns_ok, vars_)
        bad = summarize_lane_contrasts(contrasts, ns_bad, vars_)
        self.assertFalse(ok["se_undefined_lane_n_lt_2"])
        self.assertTrue(np.all(np.isfinite(ok["se"])))
        self.assertTrue(bad["se_undefined_lane_n_lt_2"])
        self.assertTrue(np.all(np.isnan(bad["se"])))
        self.assertTrue(np.all(np.isnan(bad["ci95_low"])))
        self.assertEqual(ok["significance"], SIGNIFICANCE_NOT_ESTIMABLE)
        self.assertEqual(bad["significance"], SIGNIFICANCE_NOT_ESTIMABLE)

    def test_normal_interval_uses_predeclared_1_96(self):
        contrasts = np.zeros((8, 1))
        ns = [np.full(8, 4), np.full(8, 4)]
        # each component var(mean)=1 so contrast var=2; sum=16; SE=sqrt(16/64)=0.5
        vars_ = [np.full((8, 1), 4.0), np.full((8, 1), 4.0)]
        out = summarize_lane_contrasts(contrasts, ns, vars_, z=1.96)
        np.testing.assert_allclose(out["se"], np.array([0.5]))
        np.testing.assert_allclose(out["ci95_low"], np.array([-0.98]))
        np.testing.assert_allclose(out["ci95_high"], np.array([0.98]))


class TestRidgeHoldoutNoLeakageAndIntercept(unittest.TestCase):
    def test_unpenalized_intercept_absorbs_response_shift(self):
        design = np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
        response = np.array([[1.0, 2.0], [2.0, 1.0], [3.0, 3.0]])
        weights, intercept = fit_ridge_centered(design, response, alpha=1.0)
        shifted_w, shifted_b = fit_ridge_centered(design, response + 4.0, alpha=1.0)
        np.testing.assert_allclose(shifted_w, weights, atol=1e-10)
        np.testing.assert_allclose(shifted_b, intercept + 4.0, atol=1e-10)
        pred = predict_ridge(np.zeros((1, 2)), weights, intercept)
        np.testing.assert_allclose(pred[0], intercept, atol=1e-10)

    def test_held_out_feature_and_test_outcomes_cannot_enter_fit(self):
        # Columns: A, B, canary-only-in-test
        train_x = np.array(
            [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 1.0, 0.0],
            ]
        )
        train_y = np.array(
            [
                [1.0, 0.5],
                [0.0, 1.0],
                [1.2, 1.4],
            ]
        )
        val_x = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        val_y = np.array([[0.9, 0.4], [0.1, 1.1]])
        test_x = np.array([[0.0, 0.0, 1.0]])
        test_y_small = np.array([[0.0, 0.0]])
        test_y_huge = np.array([[400.0, -400.0]])
        first = run_ridge_selection(
            train_x,
            train_y,
            val_x,
            val_y,
            alphas=(0.1, 1.0, 10.0, 100.0, 1000.0),
            test_x=test_x,
            test_y=test_y_small,
        )
        second = run_ridge_selection(
            train_x,
            train_y,
            val_x,
            val_y,
            alphas=(0.1, 1.0, 10.0, 100.0, 1000.0),
            test_x=test_x,
            test_y=test_y_huge,
        )
        self.assertEqual(first["alpha"], second["alpha"])
        np.testing.assert_allclose(first["weights"], second["weights"])
        np.testing.assert_allclose(first["intercept"], second["intercept"])
        np.testing.assert_allclose(first["train_prediction"], second["train_prediction"])
        np.testing.assert_allclose(first["weights"][2], np.array([0.0, 0.0]), atol=1e-10)
        leaked = fit_ridge_centered(
            np.vstack([train_x, test_x]),
            np.vstack([train_y, test_y_huge]),
            alpha=first["alpha"],
        )[0]
        self.assertGreater(float(np.max(np.abs(leaked[2]))), 1.0)

    def test_alpha_is_argmin_of_validation_macro_mse_only(self):
        train_x = np.array([[1.0], [1.0], [1.0], [0.0]])
        train_y = np.array([[1.0], [1.1], [0.9], [0.0]])
        val_x = np.array([[1.0], [0.0]])
        val_y = np.array([[10.0], [0.0]])
        result = run_ridge_selection(
            train_x,
            train_y,
            val_x,
            val_y,
            alphas=(0.1, 1.0, 10.0, 100.0, 1000.0),
            test_x=np.array([[1.0]]),
            test_y=np.array([[-50.0]]),
        )
        independent = []
        for alpha in (0.1, 1.0, 10.0, 100.0, 1000.0):
            weights, intercept = fit_ridge_centered(train_x, train_y, alpha)
            pred = predict_ridge(val_x, weights, intercept)
            mse = float(np.mean((pred - val_y) ** 2))
            independent.append((mse, alpha))
        expected_alpha = min(independent)[1]
        self.assertEqual(result["alpha"], expected_alpha)
        self.assertEqual(result["test_used_in_selection"], False)


class TestZeroPearsonAndCosine(unittest.TestCase):
    def test_constant_or_zero_vector_is_nan(self):
        self.assertTrue(np.isnan(delta_pearson(np.array([1.0, 1.0, 1.0]), np.array([1.0, 2.0, 3.0]))))
        self.assertTrue(np.isnan(delta_pearson(np.array([0.0, 0.0, 0.0]), np.array([0.0, 0.0, 0.0]))))
        self.assertTrue(np.isnan(cosine_similarity(np.array([0.0, 0.0, 0.0]), np.array([1.0, 2.0, 3.0]))))
        self.assertTrue(np.isnan(cosine_similarity(np.array([1.0, 2.0]), np.array([0.0, 0.0]))))

    def test_nonzero_orthogonal_cosine_is_zero_and_pearson_is_defined(self):
        self.assertAlmostEqual(cosine_similarity(np.array([1.0, 0.0]), np.array([0.0, 1.0])), 0.0)
        value = delta_pearson(np.array([0.0, 1.0, 2.0]), np.array([0.0, 2.0, 4.0]))
        self.assertFalse(np.isnan(value))
        self.assertAlmostEqual(value, 1.0)

    def test_skill_and_ratio_are_nan_when_denominator_is_zero(self):
        self.assertTrue(np.isnan(skill_score(0.2, 0.0)))
        self.assertTrue(np.isnan(zero_safe_ratio(0.3, 0.0)))
        self.assertAlmostEqual(skill_score(4.0, 8.0), 0.5)
        self.assertAlmostEqual(zero_safe_ratio(3.0, 6.0), 0.5)


class TestCanonicalOrientationGrouping(unittest.TestCase):
    def test_single_and_pair_orientations_share_one_key(self):
        labels = (
            "TSC22D1+ctrl",
            "ctrl+TSC22D1",
            "KLF1+MAP2K6",
            "MAP2K6+KLF1",
            "ctrl",
        )
        grouped = group_orientation_labels(labels)
        self.assertEqual(grouped["TSC22D1"], ("TSC22D1+ctrl", "ctrl+TSC22D1"))
        self.assertEqual(grouped["KLF1+MAP2K6"], ("KLF1+MAP2K6", "MAP2K6+KLF1"))
        self.assertEqual(grouped["ctrl"], ("ctrl",))
        self.assertEqual(canonicalize_condition("ctrl+TSC22D1"), "TSC22D1")
        self.assertEqual(len(grouped), 3)

    def test_dual_orientation_singles_are_those_with_two_labels(self):
        labels = [
            "CEBPE+ctrl",
            "ctrl+CEBPE",
            "AHR+ctrl",
            "KLF1+MAP2K6",
        ]
        grouped = group_orientation_labels(labels)
        dual = [key for key, originals in grouped.items() if "+" not in key and key != "ctrl" and len(originals) == 2]
        self.assertEqual(dual, ["CEBPE"])
        self.assertEqual(len(grouped["AHR"]), 1)


class TestTestOutcomesCannotAlterSelectedAlphaOrTrainPredictions(unittest.TestCase):
    def test_mutating_test_matrix_leaves_train_artifacts_identical(self):
        rng = np.random.default_rng(20260911)
        train_x = rng.integers(0, 2, size=(12, 5)).astype(np.float64)
        train_x[:, 0] = 1.0
        train_y = train_x @ np.array(
            [
                [1.0, 0.0, -0.5],
                [0.0, 1.0, 0.2],
                [0.3, 0.0, 0.0],
                [0.0, 0.4, 0.1],
                [0.0, 0.0, 0.0],
            ]
        ) + np.array([0.5, -0.2, 0.1])
        val_x = rng.integers(0, 2, size=(4, 5)).astype(np.float64)
        val_y = val_x @ np.array(
            [
                [1.0, 0.0, -0.5],
                [0.0, 1.0, 0.2],
                [0.3, 0.0, 0.0],
                [0.0, 0.4, 0.1],
                [0.0, 0.0, 0.0],
            ]
        )
        test_x = rng.integers(0, 2, size=(3, 5)).astype(np.float64)
        a = run_ridge_selection(
            train_x,
            train_y,
            val_x,
            val_y,
            alphas=(0.1, 1.0, 10.0, 100.0, 1000.0),
            test_x=test_x,
            test_y=np.zeros((3, 3)),
        )
        b = run_ridge_selection(
            train_x,
            train_y,
            val_x,
            val_y,
            alphas=(0.1, 1.0, 10.0, 100.0, 1000.0),
            test_x=test_x,
            test_y=rng.normal(scale=50.0, size=(3, 3)),
        )
        self.assertEqual(a["alpha"], b["alpha"])
        np.testing.assert_array_equal(a["train_prediction"], b["train_prediction"])
        self.assertEqual(a["n_train_rows"], 12)
        self.assertEqual(a["n_val_rows"], 4)

    def test_matching_single_is_average_not_sum(self):
        delta_a = np.array([2.0, 0.0])
        delta_b = np.array([0.0, 4.0])
        np.testing.assert_array_equal(matching_single_delta(delta_a, delta_b), np.array([1.0, 2.0]))
        np.testing.assert_array_equal(additive_delta(delta_a, delta_b), np.array([2.0, 4.0]))


class TestLeaveOneLaneOutAndMultiHot(unittest.TestCase):
    def test_max_rms_change_is_zero_when_lanes_identical(self):
        contrasts = np.ones((8, 3))
        summary = leave_one_lane_out_rms(contrasts)
        self.assertAlmostEqual(summary["max_rms_change"], 0.0)
        self.assertAlmostEqual(summary["full_effect_rms"], 1.0)

    def test_one_outlier_lane_changes_rms(self):
        contrasts = np.zeros((8, 1))
        contrasts[0, 0] = 8.0
        summary = leave_one_lane_out_rms(contrasts)
        self.assertGreater(summary["max_rms_change"], 0.0)
        self.assertAlmostEqual(vector_rms(np.array([1.0, 0.0, 0.0, 0.0])), 0.5)

    def test_multi_hot_encodes_main_effects_only(self):
        names = ("AHR", "KLF1", "MAP2K6")
        rows = multi_hot_rows(("AHR", "KLF1+MAP2K6", "ctrl"), names)
        np.testing.assert_array_equal(
            rows,
            np.array(
                [
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 1.0],
                    [0.0, 0.0, 0.0],
                ]
            ),
        )


class TestSparseLaneStatsDoNotNeedDenseCounts(unittest.TestCase):
    def test_csr_block_means_match_known_values(self):
        from modeling_lib import compute_group_lane_stats

        data = np.array([10.0, 0.0, 0.0, 4.0, 6.0], dtype=np.float64)
        X = sp.csr_matrix(
            np.array(
                [
                    [10.0, 1.0],
                    [0.0, 1.0],
                    [0.0, 1.0],
                    [4.0, 2.0],
                    [6.0, 2.0],
                ]
            )
        )
        stats = compute_group_lane_stats(
            X,
            group_labels=["G", "G", "G", "C", "C"],
            lanes=[1, 2, 2, 1, 1],
            lane_ids=(1, 2),
        )
        np.testing.assert_allclose(stats["G"]["mean"][0], np.array([10.0, 1.0]))
        np.testing.assert_allclose(stats["G"]["mean"][1], np.array([0.0, 1.0]))
        np.testing.assert_equal(stats["G"]["n"], np.array([1, 2]))
        self.assertTrue(np.isnan(stats["G"]["sample_var"][0, 0]))
        np.testing.assert_allclose(stats["C"]["mean"][0], np.array([5.0, 2.0]))
        self.assertEqual(data[0], 10.0)


if __name__ == "__main__":
    unittest.main()
