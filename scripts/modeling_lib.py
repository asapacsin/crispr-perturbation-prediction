"""Testable math for lane-balanced effects, additivity residuals, and ridge baselines.

No I/O. No H5AD access. Callers supply arrays. Expected values in tests are
independent literals, not copies of these formulas.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Iterable, Mapping, Sequence

import numpy as np
import scipy.sparse as sp

from dataset_schema_lib import CONTROL_TOKEN, canonicalize_condition

SIGNIFICANCE_NOT_ESTIMABLE = "not_estimable_independent_bioreplicates_absent"
LANE_IDS = (1, 2, 3, 4, 5, 6, 7, 8)
N_LANES = 8
Z_NORMAL_95 = 1.96
RIDGE_ALPHAS = (0.1, 1.0, 10.0, 100.0, 1000.0)
SPLIT_SEED = 20260911
BOOTSTRAP_SEED = 20260913
N_BOOTSTRAP = 2000


def group_orientation_labels(labels: Iterable[object]) -> dict[str, tuple[str, ...]]:
    """Map each canonical target set to its original orientation labels."""
    buckets: dict[str, list[str]] = defaultdict(list)
    seen: dict[str, set[str]] = defaultdict(set)
    for label in labels:
        key = canonicalize_condition(label)
        text = str(label)
        if text not in seen[key]:
            seen[key].add(text)
            buckets[key].append(text)
    return {key: tuple(buckets[key]) for key in sorted(buckets)}


def dual_orientation_single_keys(
    grouped: Mapping[str, Sequence[str]],
) -> list[str]:
    """Canonical single-gene keys that appear under both stored orientations."""
    out = []
    for key, originals in grouped.items():
        if key == CONTROL_TOKEN:
            continue
        if "+" in key:
            continue
        if len(originals) == 2:
            out.append(key)
    return out


def multi_hot_rows(
    condition_keys: Sequence[str],
    target_names: Sequence[str],
) -> np.ndarray:
    """Main-effect multi-hot only. Control is the zero row. No pair one-hots."""
    index = {name: i for i, name in enumerate(target_names)}
    design = np.zeros((len(condition_keys), len(target_names)), dtype=np.float64)
    for row, key in enumerate(condition_keys):
        if key == CONTROL_TOKEN:
            continue
        for token in str(key).split("+"):
            if token == CONTROL_TOKEN or token == "":
                continue
            design[row, index[token]] = 1.0
    return design


def lane_balanced_mean(lane_means: np.ndarray) -> np.ndarray:
    """Average of within-lane means. Missing lanes (NaN rows) are omitted."""
    values = np.asarray(lane_means, dtype=np.float64)
    return np.nanmean(values, axis=0)


def weighted_cell_mean(cell_values: np.ndarray) -> np.ndarray:
    """Ordinary cell-weighted mean, used only to contrast with lane balance."""
    values = np.asarray(cell_values, dtype=np.float64)
    if values.ndim == 1:
        values = values.reshape(-1, 1)
    return values.mean(axis=0)


def single_effect_lane(single_mean: np.ndarray, ctrl_mean: np.ndarray) -> np.ndarray:
    return np.asarray(single_mean, dtype=np.float64) - np.asarray(ctrl_mean, dtype=np.float64)


def pair_interaction_lane(
    pair_mean: np.ndarray,
    mean_a: np.ndarray,
    mean_b: np.ndarray,
    ctrl_mean: np.ndarray,
) -> np.ndarray:
    return (
        np.asarray(pair_mean, dtype=np.float64)
        - np.asarray(mean_a, dtype=np.float64)
        - np.asarray(mean_b, dtype=np.float64)
        + np.asarray(ctrl_mean, dtype=np.float64)
    )


def additive_delta(delta_a: np.ndarray, delta_b: np.ndarray) -> np.ndarray:
    return np.asarray(delta_a, dtype=np.float64) + np.asarray(delta_b, dtype=np.float64)


def matching_single_delta(delta_a: np.ndarray, delta_b: np.ndarray) -> np.ndarray:
    return 0.5 * (np.asarray(delta_a, dtype=np.float64) + np.asarray(delta_b, dtype=np.float64))


def lane_mean_sampling_variance(sample_var: np.ndarray, n: np.ndarray) -> np.ndarray:
    """Var(lane mean) = sample_var / n. NaN when n < 2. Do not invent variance."""
    variances = np.asarray(sample_var, dtype=np.float64)
    counts = np.asarray(n, dtype=np.float64)
    if counts.ndim == 1 and variances.ndim == 2:
        counts = counts.reshape(-1, 1)
    out = np.full(variances.shape, np.nan, dtype=np.float64)
    ok = counts >= 2
    np.divide(variances, counts, out=out, where=ok)
    return out


def se_of_averaged_lane_contrasts(lane_contrast_vars: np.ndarray) -> np.ndarray:
    """SE of the mean of L lane contrasts: sqrt(sum var / L^2). For 8 lanes, /64."""
    values = np.asarray(lane_contrast_vars, dtype=np.float64)
    n_lanes = values.shape[0]
    return np.sqrt(np.sum(values, axis=0) / (n_lanes ** 2))


def summarize_lane_contrasts(
    lane_contrasts: np.ndarray,
    component_ns: Sequence[np.ndarray],
    component_sample_vars: Sequence[np.ndarray],
    *,
    z: float = Z_NORMAL_95,
) -> dict:
    contrasts = np.asarray(lane_contrasts, dtype=np.float64)
    n_lanes = contrasts.shape[0]
    estimate = np.nanmean(contrasts, axis=0)
    invalid = np.zeros(n_lanes, dtype=bool)
    for counts in component_ns:
        invalid |= np.asarray(counts) < 2
    sampling_vars = [
        lane_mean_sampling_variance(sample_var, counts)
        for sample_var, counts in zip(component_sample_vars, component_ns)
    ]
    lane_contrast_var = np.zeros(contrasts.shape, dtype=np.float64)
    for part in sampling_vars:
        lane_contrast_var = lane_contrast_var + part
    if invalid.any():
        se = np.full(estimate.shape, np.nan, dtype=np.float64)
    else:
        se = se_of_averaged_lane_contrasts(lane_contrast_var)
    return {
        "estimate": estimate,
        "se": se,
        "ci95_low": estimate - z * se,
        "ci95_high": estimate + z * se,
        "positive_lane_fraction": np.mean(contrasts > 0, axis=0),
        "same_sign_lane_fraction": np.mean(
            np.sign(contrasts) == np.sign(estimate),
            axis=0,
        ),
        "minimum_lane_n": int(min(np.min(np.asarray(counts)) for counts in component_ns)),
        "se_undefined_lane_n_lt_2": bool(invalid.any()),
        "significance": SIGNIFICANCE_NOT_ESTIMABLE,
        "n_lanes": n_lanes,
        "n_lanes_used": int(np.sum(~np.all(np.isnan(contrasts), axis=1))),
    }


def vector_rms(values: np.ndarray) -> float:
    flat = np.asarray(values, dtype=np.float64).ravel()
    if flat.size == 0:
        return float("nan")
    return float(np.sqrt(np.mean(flat ** 2)))


def zero_safe_ratio(numerator: float, denominator: float) -> float:
    if denominator == 0 or not np.isfinite(denominator):
        return float("nan")
    return float(numerator) / float(denominator)


def skill_score(model_loss: float, baseline_loss: float) -> float:
    if baseline_loss == 0 or not np.isfinite(baseline_loss):
        return float("nan")
    return 1.0 - float(model_loss) / float(baseline_loss)


def _finite_std_is_zero(values: np.ndarray) -> bool:
    return float(np.std(values, ddof=0)) == 0.0


def delta_pearson(left: np.ndarray, right: np.ndarray) -> float:
    a = np.asarray(left, dtype=np.float64).ravel()
    b = np.asarray(right, dtype=np.float64).ravel()
    if a.size == 0 or b.size == 0 or a.size != b.size:
        return float("nan")
    if _finite_std_is_zero(a) or _finite_std_is_zero(b):
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def cosine_similarity(left: np.ndarray, right: np.ndarray) -> float:
    a = np.asarray(left, dtype=np.float64).ravel()
    b = np.asarray(right, dtype=np.float64).ravel()
    norm_a = float(np.linalg.norm(a))
    norm_b = float(np.linalg.norm(b))
    if norm_a == 0.0 or norm_b == 0.0:
        return float("nan")
    return float(np.dot(a, b) / (norm_a * norm_b))


def per_condition_mse(pred: np.ndarray, truth: np.ndarray) -> np.ndarray:
    residual = np.asarray(pred, dtype=np.float64) - np.asarray(truth, dtype=np.float64)
    return np.mean(residual ** 2, axis=1)


def per_condition_mae(pred: np.ndarray, truth: np.ndarray) -> np.ndarray:
    residual = np.asarray(pred, dtype=np.float64) - np.asarray(truth, dtype=np.float64)
    return np.mean(np.abs(residual), axis=1)


def macro_mse(pred: np.ndarray, truth: np.ndarray) -> float:
    return float(np.mean(per_condition_mse(pred, truth)))


def macro_mae(pred: np.ndarray, truth: np.ndarray) -> float:
    return float(np.mean(per_condition_mae(pred, truth)))


def leave_one_lane_out_rms(lane_contrasts: np.ndarray) -> dict:
    contrasts = np.asarray(lane_contrasts, dtype=np.float64)
    n_lanes = contrasts.shape[0]
    full = np.nanmean(contrasts, axis=0)
    lolo_rms = np.empty(n_lanes, dtype=np.float64)
    rms_change = np.empty(n_lanes, dtype=np.float64)
    for left_out in range(n_lanes):
        others = np.delete(contrasts, left_out, axis=0)
        estimate = np.nanmean(others, axis=0)
        lolo_rms[left_out] = vector_rms(estimate)
        rms_change[left_out] = vector_rms(estimate - full)
    return {
        "full": full,
        "max_rms_change": float(np.max(rms_change)),
        "min_lolo_effect_rms": float(np.min(lolo_rms)),
        "max_lolo_effect_rms": float(np.max(lolo_rms)),
        "full_effect_rms": vector_rms(full),
    }


def compute_group_lane_stats(
    X,
    group_labels: Sequence[str],
    lanes: Sequence[int],
    lane_ids: Sequence[int] = LANE_IDS,
) -> dict[str, dict[str, np.ndarray]]:
    """Per-group, per-lane cell count, mean, and sample variance on sparse X."""
    matrix = X.tocsr() if sp.issparse(X) else sp.csr_matrix(np.asarray(X))
    n_genes = matrix.shape[1]
    lane_index = {int(lane): i for i, lane in enumerate(lane_ids)}
    n_lanes = len(lane_ids)
    group_arr = np.asarray(group_labels)
    lane_arr = np.asarray(lanes, dtype=int)
    result: dict[str, dict[str, np.ndarray]] = {}
    for group in sorted(set(group_arr.tolist()), key=str):
        counts = np.zeros(n_lanes, dtype=np.int64)
        means = np.full((n_lanes, n_genes), np.nan, dtype=np.float64)
        variances = np.full((n_lanes, n_genes), np.nan, dtype=np.float64)
        in_group = group_arr == group
        for lane, slot in lane_index.items():
            rows = np.flatnonzero(in_group & (lane_arr == lane))
            counts[slot] = rows.size
            if rows.size == 0:
                continue
            block = matrix[rows].toarray().astype(np.float64, copy=False)
            means[slot] = block.mean(axis=0)
            if rows.size >= 2:
                variances[slot] = block.var(axis=0, ddof=1)
        result[str(group)] = {"n": counts, "mean": means, "sample_var": variances}
    return result


def fit_ridge_centered(
    design: np.ndarray,
    response: np.ndarray,
    alpha: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Ridge with unpenalized intercept via column-centering design and response."""
    x = np.asarray(design, dtype=np.float64)
    y = np.asarray(response, dtype=np.float64)
    if y.ndim == 1:
        y = y.reshape(-1, 1)
    x_mean = x.mean(axis=0)
    y_mean = y.mean(axis=0)
    x_c = x - x_mean
    y_c = y - y_mean
    gram = x_c.T @ x_c
    n_features = gram.shape[0]
    gram.flat[:: n_features + 1] += float(alpha)
    weights = np.linalg.solve(gram, x_c.T @ y_c)
    intercept = y_mean - x_mean @ weights
    return weights, intercept


def predict_ridge(design: np.ndarray, weights: np.ndarray, intercept: np.ndarray) -> np.ndarray:
    x = np.asarray(design, dtype=np.float64)
    return x @ np.asarray(weights, dtype=np.float64) + np.asarray(intercept, dtype=np.float64)


def run_ridge_selection(
    train_x: np.ndarray,
    train_y: np.ndarray,
    val_x: np.ndarray,
    val_y: np.ndarray,
    alphas: Sequence[float] = RIDGE_ALPHAS,
    test_x: np.ndarray | None = None,
    test_y: np.ndarray | None = None,
) -> dict:
    """Select one alpha on validation macro MSE. Test outcomes are ignored."""
    del test_y  # protocol: test cannot alter selection or the train fit
    curve = []
    best: dict | None = None
    for alpha in alphas:
        weights, intercept = fit_ridge_centered(train_x, train_y, alpha)
        val_pred = predict_ridge(val_x, weights, intercept)
        mse = macro_mse(val_pred, val_y)
        curve.append({"alpha": float(alpha), "validation_macro_mse": mse})
        better = best is None or mse < best["validation_macro_mse"] or (
            mse == best["validation_macro_mse"] and float(alpha) < float(best["alpha"])
        )
        if better:
            best = {
                "alpha": float(alpha),
                "validation_macro_mse": mse,
                "weights": weights,
                "intercept": intercept,
            }
    assert best is not None
    train_prediction = predict_ridge(train_x, best["weights"], best["intercept"])
    payload = {
        "alpha": best["alpha"],
        "weights": best["weights"],
        "intercept": best["intercept"],
        "train_prediction": train_prediction,
        "validation_prediction": predict_ridge(val_x, best["weights"], best["intercept"]),
        "validation_curve": curve,
        "n_train_rows": int(np.asarray(train_x).shape[0]),
        "n_val_rows": int(np.asarray(val_x).shape[0]),
        "test_used_in_selection": False,
    }
    if test_x is not None:
        payload["test_prediction"] = predict_ridge(test_x, best["weights"], best["intercept"])
    return payload


def paired_bootstrap_mean_difference(
    model_scores: np.ndarray,
    baseline_scores: np.ndarray,
    *,
    seed: int = BOOTSTRAP_SEED,
    n_samples: int = N_BOOTSTRAP,
) -> dict:
    """Condition-level paired bootstrap. Conditional on a fixed model and reference."""
    model = np.asarray(model_scores, dtype=np.float64).ravel()
    baseline = np.asarray(baseline_scores, dtype=np.float64).ravel()
    if model.size != baseline.size:
        raise ValueError("bootstrap score vectors must align")
    diffs = model - baseline
    rng = np.random.default_rng(seed)
    n = diffs.size
    samples = np.empty(n_samples, dtype=np.float64)
    for i in range(n_samples):
        take = rng.integers(0, n, n)
        samples[i] = diffs[take].mean()
    return {
        "mean_difference": float(diffs.mean()),
        "bootstrap_mean": float(samples.mean()),
        "bootstrap_std": float(samples.std(ddof=1)),
        "percentile_2_5": float(np.percentile(samples, 2.5)),
        "percentile_97_5": float(np.percentile(samples, 97.5)),
        "n_samples": int(n_samples),
        "seed": int(seed),
        "n_conditions": int(n),
        "interpretation": "conditional_on_fixed_model_and_train_reference",
        "not_biological_confidence_interval": True,
    }


def condition_signature_svd(matrix: np.ndarray, n_components: int = 10) -> dict:
    """Column-centered SVD of condition signatures. Descriptive, not a training feature."""
    values = np.asarray(matrix, dtype=np.float64)
    centered = values - values.mean(axis=0, keepdims=True)
    u, singular, vt = np.linalg.svd(centered, full_matrices=False)
    k = min(n_components, singular.size)
    return {
        "scores": u[:, :k] * singular[:k],
        "loadings": vt[:k],
        "singular_values": singular[:k],
        "n_components": k,
    }


def pair_partners(canonical_key: str) -> tuple[str, str]:
    tokens = [token for token in canonical_key.split("+") if token and token != CONTROL_TOKEN]
    if len(tokens) != 2:
        raise ValueError(f"expected a two-gene canonical key, got {canonical_key!r}")
    return tokens[0], tokens[1]


def json_safe(value):
    if isinstance(value, (np.floating, float)):
        number = float(value)
        if not np.isfinite(number):
            return None
        return number
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, np.ndarray):
        return [json_safe(item) for item in value.tolist()]
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    return value
