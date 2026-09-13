"""Descriptive diagnostics on frozen Package-B outputs.

Reads existing lane_stats.h5, genes.csv, describe tables, and predictions.npz.
Does not reload the H5AD, refit, retune, or touch the main modeling report.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

try:
    import yaml
except ImportError as exc:  # pragma: no cover - environment check
    raise SystemExit(
        "ERROR: PyYAML is required to read config/dataset_schema.yaml. "
        "This script will not pip-install it. "
        f"Original import error: {exc}"
    ) from exc

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

from modeling_lib import (
    SIGNIFICANCE_NOT_ESTIMABLE,
    Z_NORMAL_95,
    cosine_similarity,
    delta_pearson,
    json_safe,
    lane_balanced_mean,
    macro_mae,
    macro_mse,
    pair_interaction_lane,
    pair_partners,
    single_effect_lane,
    summarize_lane_contrasts,
    vector_rms,
    zero_safe_ratio,
)

MODELING_DIR = ROOT / "reports" / "modeling"
PREPARE_DIR = MODELING_DIR / "prepare"
FIT_DIR = MODELING_DIR / "fit"
DESCRIBE_DIR = MODELING_DIR / "describe"
DIAG_DIR = MODELING_DIR / "diagnostics"
LANE_STATS_PATH = PREPARE_DIR / "lane_stats.h5"
REPORT_PATH = ROOT / "reports" / "MODELING_DIAGNOSTICS.md"
RESOURCE_EVIDENCE = ROOT / "reports" / "modeling_resource_evidence.json"
MIN_LANE_SIGN = 7.0 / 8.0
REQUIRED = (
    LANE_STATS_PATH,
    PREPARE_DIR / "genes.csv",
    PREPARE_DIR / "condition_sizes.csv",
    PREPARE_DIR / "prepare_meta.json",
    FIT_DIR / "predictions.npz",
    FIT_DIR / "summary_metrics.json",
    FIT_DIR / "per_condition_metrics.csv",
    FIT_DIR / "train_conditions.csv",
    DESCRIBE_DIR / "pair_summaries.csv",
    DESCRIBE_DIR / "svd_modules.csv",
)


class DiagnosticsError(RuntimeError):
    pass


def map_own_target_indices(
    target_names: tuple[str, ...] | list[str],
    gene_names: tuple[str, ...] | list[str] | np.ndarray,
) -> tuple[dict[str, int | None], list[str]]:
    """Map stored target symbols onto unique measured gene_name indices."""
    names = [str(name) for name in gene_names]
    counts: dict[str, list[int]] = {}
    for index, name in enumerate(names):
        counts.setdefault(name, []).append(index)
    mapping: dict[str, int | None] = {}
    missing: list[str] = []
    for target in target_names:
        hits = counts.get(str(target), [])
        if len(hits) != 1:
            mapping[str(target)] = None
            missing.append(str(target))
        else:
            mapping[str(target)] = hits[0]
    return mapping, missing


def own_target_value(effect_vector: np.ndarray, index: int | None) -> float:
    if index is None:
        return float("nan")
    return float(np.asarray(effect_vector, dtype=np.float64)[int(index)])


def project_residual_onto_shift(residual: np.ndarray, direction: np.ndarray) -> dict[str, float]:
    residual_v = np.asarray(residual, dtype=np.float64).ravel()
    shift = np.asarray(direction, dtype=np.float64).ravel()
    shift_ssq = float(np.dot(shift, shift))
    resid_ssq = float(np.dot(residual_v, residual_v))
    nan = {
        "coefficient": float("nan"),
        "squared_norm_fraction": float("nan"),
        "cosine_to_neg_direction": float("nan"),
        "orthogonal_rms": float("nan"),
    }
    if shift_ssq == 0.0 or resid_ssq == 0.0:
        return nan
    coefficient = float(np.dot(residual_v, shift) / shift_ssq)
    projection = coefficient * shift
    fraction = float(np.dot(projection, projection) / resid_ssq)
    orthogonal = residual_v - projection
    return {
        "coefficient": coefficient,
        "squared_norm_fraction": fraction,
        "cosine_to_neg_direction": cosine_similarity(residual_v, -shift),
        "orthogonal_rms": vector_rms(orthogonal),
    }


def top_cosine_pairs(
    signatures: np.ndarray,
    labels: list[str] | tuple[str, ...],
    n_pairs: int = 10,
) -> list[dict]:
    matrix = np.asarray(signatures, dtype=np.float64)
    names = [str(label) for label in labels]
    n = matrix.shape[0]
    scored: list[tuple[float, int, int]] = []
    for i in range(n):
        for j in range(i + 1, n):
            scored.append((cosine_similarity(matrix[i], matrix[j]), i, j))
    scored.sort(key=lambda row: (-row[0] if np.isfinite(row[0]) else 1.0, names[row[1]], names[row[2]]))
    out = []
    for rank, (value, i, j) in enumerate(scored[:n_pairs], start=1):
        out.append(
            {
                "rank": rank,
                "gene_a": names[i],
                "gene_b": names[j],
                "pair": (names[i], names[j]),
                "cosine": value,
            }
        )
    return out


def select_example_pairs(
    pair_keys: list[str],
    residual_rms: np.ndarray,
    observed_rms: np.ndarray,
    relative_residual: np.ndarray,
    top3_min_sign: np.ndarray,
    min_sign: float = MIN_LANE_SIGN,
) -> list[dict]:
    keys = list(pair_keys)
    resid = np.asarray(residual_rms, dtype=np.float64)
    observed = np.asarray(observed_rms, dtype=np.float64)
    relative = np.asarray(relative_residual, dtype=np.float64)
    signs = np.asarray(top3_min_sign, dtype=np.float64)
    chosen: list[dict] = []
    used: set[str] = set()

    largest = int(np.nanargmax(resid))
    chosen.append(
        {
            "canonical_perturbation": keys[largest],
            "role": "largest_canonical_residual",
            "index": largest,
        }
    )
    used.add(keys[largest])

    median_obs = float(np.nanmedian(observed))
    eligible = [
        i
        for i, key in enumerate(keys)
        if key not in used and observed[i] > median_obs and np.isfinite(relative[i])
    ]
    if eligible:
        best = min(eligible, key=lambda i: (relative[i], -observed[i], keys[i]))
        chosen.append(
            {
                "canonical_perturbation": keys[best],
                "role": "lowest_relative_residual_above_median_observed_rms",
                "index": best,
                "median_observed_rms": median_obs,
            }
        )
        used.add(keys[best])

    signed = [
        i
        for i, key in enumerate(keys)
        if key not in used and np.isfinite(signs[i]) and signs[i] >= min_sign
    ]
    if signed:
        best = max(signed, key=lambda i: (resid[i], keys[i]))
        chosen.append(
            {
                "canonical_perturbation": keys[best],
                "role": "largest_lane_sign_supported_residual",
                "index": best,
                "min_same_sign_top3": float(signs[best]),
            }
        )
    return chosen


def read_lane_stats_group(container: h5py.Group) -> dict:
    out = {}
    for key in container.keys():
        out[key] = {
            "n": container[key]["n"][()],
            "mean": container[key]["mean"][()],
            "sample_var": container[key]["sample_var"][()],
        }
    return out


def write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(json_safe(payload), indent=2, sort_keys=True),
        encoding="utf-8",
    )


def _fmt(value) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        if not np.isfinite(value):
            return "NA"
        return f"{value:.6g}"
    return str(value)


def _md_table(rows: list[dict], columns: list[str]) -> str:
    header = "| " + " | ".join(columns) + " |"
    sep = "| " + " | ".join("---" for _ in columns) + " |"
    body = [
        "| " + " | ".join(_fmt(row.get(col, "")) for col in columns) + " |"
        for row in rows
    ]
    return "\n".join([header, sep, *body])


def require_inputs() -> None:
    missing = [str(path.relative_to(ROOT)).replace("\\", "/") for path in REQUIRED if not path.is_file()]
    if missing:
        raise DiagnosticsError(
            "required frozen artifacts are not ready: " + ", ".join(missing)
        )


def condition_mean(entry: dict) -> np.ndarray:
    return lane_balanced_mean(entry["mean"])


def recompute_split_metrics(pred: np.ndarray, truth: np.ndarray) -> dict:
    residual = np.asarray(pred, dtype=np.float64) - np.asarray(truth, dtype=np.float64)
    mse = np.mean(residual ** 2, axis=1)
    mae = np.mean(np.abs(residual), axis=1)
    pears = np.array([delta_pearson(pred[i], truth[i]) for i in range(len(pred))], dtype=np.float64)
    coss = np.array([cosine_similarity(pred[i], truth[i]) for i in range(len(pred))], dtype=np.float64)
    return {
        "macro_mse": macro_mse(pred, truth),
        "macro_mae": macro_mae(pred, truth),
        "macro_delta_pearson": float(np.nanmean(pears)) if np.any(np.isfinite(pears)) else None,
        "macro_cosine": float(np.nanmean(coss)) if np.any(np.isfinite(coss)) else None,
        "n_conditions": int(len(mse)),
        "n_pearson_defined": int(np.isfinite(pears).sum()),
        "n_cosine_defined": int(np.isfinite(coss).sum()),
        "per_condition_mse": mse,
        "per_condition_mae": mae,
    }


def _close(saved, recomputed, atol: float = 1e-10, rtol: float = 1e-8) -> bool:
    if saved is None and recomputed is None:
        return True
    if saved is None or recomputed is None:
        return False
    return bool(np.isclose(float(saved), float(recomputed), atol=atol, rtol=rtol))


def run_diagnostics() -> dict:
    require_inputs()
    DIAG_DIR.mkdir(parents=True, exist_ok=True)
    schema = yaml.safe_load((ROOT / "config" / "dataset_schema.yaml").read_text(encoding="utf-8"))
    target_names = list(schema["targets"]["names"])
    genes = pd.read_csv(PREPARE_DIR / "genes.csv")
    sizes = pd.read_csv(PREPARE_DIR / "condition_sizes.csv")
    prepare_meta = json.loads((PREPARE_DIR / "prepare_meta.json").read_text(encoding="utf-8"))
    pair_summary = pd.read_csv(DESCRIBE_DIR / "pair_summaries.csv")
    svd_modules = pd.read_csv(DESCRIBE_DIR / "svd_modules.csv")
    saved_metrics = json.loads((FIT_DIR / "summary_metrics.json").read_text(encoding="utf-8"))
    train_cond_csv = pd.read_csv(FIT_DIR / "train_conditions.csv")
    go_evidence = {}
    if RESOURCE_EVIDENCE.is_file():
        go_evidence = json.loads(RESOURCE_EVIDENCE.read_text(encoding="utf-8"))

    originals = {
        str(row.canonical_perturbation): [
            part for part in str(row.original_condition_labels).split("|") if part
        ]
        for row in sizes.itertuples(index=False)
    }
    mapping, missing_targets = map_own_target_indices(target_names, genes["gene_name"].tolist())

    with h5py.File(LANE_STATS_PATH, "r") as store:
        canonical = read_lane_stats_group(store["canonical"])
        original = read_lane_stats_group(store["original"])
    ctrl = canonical["ctrl"]
    ctrl_mean = condition_mean(ctrl)

    predictions = np.load(FIT_DIR / "predictions.npz", allow_pickle=True)
    mu0 = np.asarray(predictions["mu0"], dtype=np.float64)
    train_keys = [str(key) for key in predictions["train_keys"].tolist()]
    gene_ids = [str(item) for item in predictions["gene_ids"].tolist()]

    if gene_ids != genes["gene_id"].astype(str).tolist():
        raise DiagnosticsError("predictions.npz gene_id order does not match genes.csv")
    if train_keys != train_cond_csv["canonical_perturbation"].astype(str).tolist():
        raise DiagnosticsError("predictions.npz train_keys do not match train_conditions.csv")

    single_keys = [key for key in target_names]
    if any(key not in canonical for key in single_keys):
        raise DiagnosticsError("a schema target is missing from canonical lane_stats")

    single_effects = {}
    single_pack = {}
    own_rows = []
    train_single_deltas = []
    for key in single_keys:
        entry = canonical[key]
        contrasts = single_effect_lane(entry["mean"], ctrl["mean"])
        packed = summarize_lane_contrasts(
            contrasts,
            [entry["n"], ctrl["n"]],
            [entry["sample_var"], ctrl["sample_var"]],
            z=Z_NORMAL_95,
        )
        single_effects[key] = packed["estimate"]
        single_pack[key] = packed
        train_single_deltas.append(condition_mean(entry) - mu0)
        index = mapping[key]
        effect = own_target_value(packed["estimate"], index)
        se = float("nan") if index is None else float(packed["se"][index])
        own_rows.append(
            {
                "target_gene": key,
                "gene_id": None if index is None else str(genes.loc[index, "gene_id"]),
                "gene_index": index,
                "effect": effect,
                "conditional_cell_sampling_se": se,
                "ci95_low": None if index is None else float(packed["ci95_low"][index]),
                "ci95_high": None if index is None else float(packed["ci95_high"][index]),
                "positive_lane_fraction": None if index is None else float(packed["positive_lane_fraction"][index]),
                "same_sign_lane_fraction": None if index is None else float(packed["same_sign_lane_fraction"][index]),
                "sign": (
                    "missing"
                    if index is None
                    else "positive"
                    if effect > 0
                    else "negative"
                    if effect < 0
                    else "zero"
                ),
                "minimum_lane_n": packed["minimum_lane_n"],
                "se_undefined_lane_n_lt_2": packed["se_undefined_lane_n_lt_2"],
                "significance": SIGNIFICANCE_NOT_ESTIMABLE,
                "interpretation": "normalized_own_target_effect_not_molecular_validation",
            }
        )
    own_frame = pd.DataFrame(own_rows)
    own_counts = {
        "n_targets": len(single_keys),
        "n_mapped": int((own_frame["sign"] != "missing").sum()),
        "n_missing": len(missing_targets),
        "missing_targets": missing_targets,
        "n_positive": int((own_frame["sign"] == "positive").sum()),
        "n_negative": int((own_frame["sign"] == "negative").sum()),
        "n_zero": int((own_frame["sign"] == "zero").sum()),
        "interpretation": "counts_of_normalized_own_target_log1p_effect_not_guide_validation",
    }

    shared_shift = np.mean(np.stack(train_single_deltas), axis=0)
    pair_keys = [key for key in sorted(canonical) if "+" in key]
    if len(pair_keys) != 131:
        raise DiagnosticsError(f"expected 131 pairs in lane_stats, got {len(pair_keys)}")

    pair_obs = {}
    pair_add = {}
    pair_resid = {}
    pair_pack = {}
    projection_rows = []
    orientation_pair_rows = []
    top3_min_sign = []
    residual_rms = []
    observed_rms = []
    relative_residual = []

    dual = set(prepare_meta.get("dual_orientation_singles", []))
    for key in pair_keys:
        gene_a, gene_b = pair_partners(key)
        entry = canonical[key]
        contrasts = pair_interaction_lane(
            entry["mean"],
            canonical[gene_a]["mean"],
            canonical[gene_b]["mean"],
            ctrl["mean"],
        )
        packed = summarize_lane_contrasts(
            contrasts,
            [entry["n"], canonical[gene_a]["n"], canonical[gene_b]["n"], ctrl["n"]],
            [
                entry["sample_var"],
                canonical[gene_a]["sample_var"],
                canonical[gene_b]["sample_var"],
                ctrl["sample_var"],
            ],
            z=Z_NORMAL_95,
        )
        observed = condition_mean(entry) - ctrl_mean
        additive = single_effects[gene_a] + single_effects[gene_b]
        residual = packed["estimate"]
        pair_obs[key] = observed
        pair_add[key] = additive
        pair_resid[key] = residual
        pair_pack[key] = packed
        proj = project_residual_onto_shift(residual, shared_shift)
        obs_rms = vector_rms(observed)
        resid_rms = vector_rms(residual)
        projection_rows.append(
            {
                "canonical_perturbation": key,
                "target_a": gene_a,
                "target_b": gene_b,
                "observed_rms": obs_rms,
                "residual_rms": resid_rms,
                "residual_over_response_rms": zero_safe_ratio(resid_rms, obs_rms),
                "projection_coefficient": proj["coefficient"],
                "squared_norm_fraction": proj["squared_norm_fraction"],
                "cosine_to_neg_shared_shift": proj["cosine_to_neg_direction"],
                "orthogonal_residual_rms": proj["orthogonal_rms"],
                "role": "descriptive_alignment_not_confounding_proof",
            }
        )
        residual_rms.append(resid_rms)
        observed_rms.append(obs_rms)
        relative_residual.append(zero_safe_ratio(resid_rms, obs_rms))
        top3 = np.argsort(-np.abs(residual))[:3]
        top3_min_sign.append(float(np.min(packed["same_sign_lane_fraction"][top3])))

        labels_a = originals.get(gene_a, [gene_a])
        labels_b = originals.get(gene_b, [gene_b])
        if gene_a in dual or gene_b in dual:
            combo_rms = []
            for label_a in labels_a:
                for label_b in labels_b:
                    if label_a not in original or label_b not in original:
                        continue
                    effect_a = lane_balanced_mean(
                        single_effect_lane(original[label_a]["mean"], ctrl["mean"])
                    )
                    effect_b = lane_balanced_mean(
                        single_effect_lane(original[label_b]["mean"], ctrl["mean"])
                    )
                    alt_rms = vector_rms(observed - (effect_a + effect_b))
                    combo_rms.append(alt_rms)
                    orientation_pair_rows.append(
                        {
                            "canonical_perturbation": key,
                            "label_a": label_a,
                            "label_b": label_b,
                            "residual_rms": alt_rms,
                            "canonical_residual_rms": resid_rms,
                            "n_combinations": len(labels_a) * len(labels_b),
                        }
                    )
            if combo_rms:
                orientation_pair_rows.append(
                    {
                        "canonical_perturbation": key,
                        "label_a": "RANGE",
                        "label_b": "RANGE",
                        "residual_rms_min": min(combo_rms),
                        "residual_rms_max": max(combo_rms),
                        "canonical_residual_rms": resid_rms,
                        "n_combinations": len(combo_rms),
                        "do_not_choose_favorable_orientation": True,
                    }
                )

    projection = pd.DataFrame(projection_rows).sort_values(
        "squared_norm_fraction",
        ascending=False,
    )
    residual_rms_arr = np.asarray(residual_rms, dtype=np.float64)
    observed_rms_arr = np.asarray(observed_rms, dtype=np.float64)
    relative_arr = np.asarray(relative_residual, dtype=np.float64)
    sign_arr = np.asarray(top3_min_sign, dtype=np.float64)
    shift_summary = {
        "n_train_singles": 105,
        "mu0_source": "predictions.npz",
        "residual_source": "all_data_lane_stats_ctrl",
        "median_squared_norm_fraction": float(np.nanmedian(projection["squared_norm_fraction"])),
        "median_cosine_to_neg_shared_shift": float(np.nanmedian(projection["cosine_to_neg_shared_shift"])),
        "median_orthogonal_residual_rms": float(np.nanmedian(projection["orthogonal_residual_rms"])),
        "largest_squared_norm_fraction_pair": str(
            projection.iloc[0]["canonical_perturbation"]
        ),
        "largest_squared_norm_fraction": float(projection.iloc[0]["squared_norm_fraction"]),
        "interpretation": (
            "High squared-norm fraction is alignment of a descriptive residual with the "
            "mean training single-delta direction. It is not proof of confounding and "
            "is not a corrected prediction or a new test metric."
        ),
    }

    orientation_single_rows = []
    for key in sorted(dual):
        labels = originals.get(key, [])
        if len(labels) != 2:
            continue
        if any(label not in original for label in labels):
            continue
        mean_a = condition_mean(original[labels[0]])
        mean_b = condition_mean(original[labels[1]])
        profile_rms = vector_rms(mean_a - mean_b)
        pooled_rms = vector_rms(single_effects[key])
        orientation_single_rows.append(
            {
                "canonical_perturbation": key,
                "label_a": labels[0],
                "label_b": labels[1],
                "profile_rms": profile_rms,
                "pooled_single_effect_rms": pooled_rms,
                "profile_over_pooled_single_rms": zero_safe_ratio(profile_rms, pooled_rms),
                "n_cells_a": int(original[labels[0]]["n"].sum()),
                "n_cells_b": int(original[labels[1]]["n"].sum()),
                "no_new_biology_claim": True,
            }
        )
    orientation_singles = pd.DataFrame(orientation_single_rows)
    orientation_pairs = pd.DataFrame(orientation_pair_rows)

    ranking_rows = []
    order_resid = np.argsort(-residual_rms_arr)
    for rank, index in enumerate(order_resid[:10], start=1):
        ranking_rows.append(
            {
                "list": "largest_residual_rms",
                "rank": rank,
                "canonical_perturbation": pair_keys[index],
                "residual_rms": residual_rms_arr[index],
                "observed_rms": observed_rms_arr[index],
                "residual_over_response_rms": relative_arr[index],
            }
        )
    finite_rel = [(relative_arr[i], -observed_rms_arr[i], i) for i in range(len(pair_keys)) if np.isfinite(relative_arr[i])]
    finite_rel.sort()
    for rank, (_, _, index) in enumerate(finite_rel[:10], start=1):
        ranking_rows.append(
            {
                "list": "lowest_residual_over_response_rms",
                "rank": rank,
                "canonical_perturbation": pair_keys[index],
                "residual_rms": residual_rms_arr[index],
                "observed_rms": observed_rms_arr[index],
                "residual_over_response_rms": relative_arr[index],
                "note": "weak observed_rms must not be called additive",
            }
        )
    rankings = pd.DataFrame(ranking_rows)

    example_slots = select_example_pairs(
        pair_keys,
        residual_rms_arr,
        observed_rms_arr,
        relative_arr,
        sign_arr,
    )
    example_gene_rows = []
    example_meta = []
    for slot in example_slots:
        key = slot["canonical_perturbation"]
        gene_a, gene_b = pair_partners(key)
        residual = pair_resid[key]
        packed = pair_pack[key]
        top5 = np.argsort(-np.abs(residual))[:5]
        example_meta.append(
            {
                **slot,
                "residual_rms": float(residual_rms_arr[slot["index"]]),
                "observed_rms": float(observed_rms_arr[slot["index"]]),
                "residual_over_response_rms": float(relative_arr[slot["index"]]),
                "top3_min_same_sign": float(sign_arr[slot["index"]]),
            }
        )
        for rank, index in enumerate(top5, start=1):
            example_gene_rows.append(
                {
                    "canonical_perturbation": key,
                    "role": slot["role"],
                    "rank": rank,
                    "gene_id": str(genes.loc[index, "gene_id"]),
                    "gene_name": str(genes.loc[index, "gene_name"]),
                    "effect_a": float(single_effects[gene_a][index]),
                    "effect_b": float(single_effects[gene_b][index]),
                    "observed_ab": float(pair_obs[key][index]),
                    "additive": float(pair_add[key][index]),
                    "residual": float(residual[index]),
                    "conditional_cell_sampling_se": float(packed["se"][index]),
                    "ci95_low": float(packed["ci95_low"][index]),
                    "ci95_high": float(packed["ci95_high"][index]),
                    "same_sign_lane_fraction": float(packed["same_sign_lane_fraction"][index]),
                    "significance": SIGNIFICANCE_NOT_ESTIMABLE,
                }
            )
    examples = pd.DataFrame(example_gene_rows)

    signature_matrix = np.stack([single_effects[key] for key in single_keys])
    neighbors = pd.DataFrame(top_cosine_pairs(signature_matrix, single_keys, n_pairs=10))
    if "pair" in neighbors.columns:
        neighbors = neighbors.drop(columns=["pair"])

    svd_signed_rows = []
    for component, block in svd_modules.groupby("component", sort=True):
        pos = block.loc[block["loading"] > 0].sort_values("loading", ascending=False)
        neg = block.loc[block["loading"] < 0].sort_values("loading", ascending=True)
        svd_signed_rows.append(
            {
                "component": int(component),
                "pole": "positive",
                "genes": "|".join(pos["gene_name"].astype(str)),
                "loadings": "|".join(f"{float(v):.6g}" for v in pos["loading"]),
                "singular_value": float(block["singular_value"].iloc[0]),
                "interpretation": "data_only_signed_loading_not_named_pathway",
            }
        )
        svd_signed_rows.append(
            {
                "component": int(component),
                "pole": "negative",
                "genes": "|".join(neg["gene_name"].astype(str)),
                "loadings": "|".join(f"{float(v):.6g}" for v in neg["loading"]),
                "singular_value": float(block["singular_value"].iloc[0]),
                "interpretation": "data_only_signed_loading_not_named_pathway",
            }
        )
    svd_signed = pd.DataFrame(svd_signed_rows)

    metric_checks = verify_saved_metrics(predictions, saved_metrics, train_keys)
    go_skip = {
        "named_go_enrichment": "skipped",
        "reason": "local_go_graph_has_no_named_terms",
        "source": str(RESOURCE_EVIDENCE.relative_to(ROOT)).replace("\\", "/")
        if RESOURCE_EVIDENCE.is_file()
        else None,
        "edges": go_evidence.get("rows"),
        "nodes": go_evidence.get("nodes"),
        "columns": go_evidence.get("columns"),
        "has_named_gene_sets": go_evidence.get("has_named_gene_sets", False),
    }

    own_frame.to_csv(DIAG_DIR / "own_target_single_effects.csv", index=False, encoding="utf-8")
    write_json(DIAG_DIR / "own_target_counts.json", own_counts)
    projection.to_csv(DIAG_DIR / "shared_shift_projections.csv", index=False, encoding="utf-8")
    write_json(DIAG_DIR / "shared_shift_summary.json", shift_summary)
    orientation_singles.to_csv(
        DIAG_DIR / "orientation_profile_rms.csv",
        index=False,
        encoding="utf-8",
    )
    orientation_pairs.to_csv(
        DIAG_DIR / "orientation_pair_additive_range.csv",
        index=False,
        encoding="utf-8",
    )
    rankings.to_csv(DIAG_DIR / "pair_rankings.csv", index=False, encoding="utf-8")
    neighbors.to_csv(DIAG_DIR / "nearest_neighbor_singles.csv", index=False, encoding="utf-8")
    svd_signed.to_csv(DIAG_DIR / "svd_modules_signed.csv", index=False, encoding="utf-8")
    examples.to_csv(DIAG_DIR / "example_pair_genes.csv", index=False, encoding="utf-8")
    write_json(DIAG_DIR / "example_pairs.json", example_meta)
    write_json(DIAG_DIR / "metric_verification.json", metric_checks)
    write_json(DIAG_DIR / "go_skip.json", go_skip)
    write_json(
        DIAG_DIR / "diagnostics_meta.json",
        {
            "refit_attempted": False,
            "h5ad_reloaded": False,
            "test_retuned": False,
            "plots_added": False,
            "named_go_enrichment": "skipped",
            "test_pairs_consumed": True,
            "next_model_requires_nested_or_external_holdout": True,
            "n_example_pairs": len(example_meta),
            "pair_summary_file_used_for_display_only": True,
            "n_pair_summary_rows": int(len(pair_summary)),
        },
    )
    write_diagnostics_report(
        own_counts,
        own_frame,
        shift_summary,
        projection,
        orientation_singles,
        orientation_pairs,
        rankings,
        example_meta,
        examples,
        neighbors,
        svd_signed,
        metric_checks,
        go_skip,
        pair_summary,
    )
    return {"own_counts": own_counts, "shift_summary": shift_summary, "verification": metric_checks}


def verify_saved_metrics(predictions, saved_metrics: dict, train_keys: list[str]) -> dict:
    train_true = np.asarray(predictions["train_true"], dtype=np.float64)
    train_pred = np.asarray(predictions["train_pred_ridge"], dtype=np.float64)
    checks = {
        "refit_attempted": False,
        "train_ridge_shape": list(train_pred.shape),
        "train_true_shape": list(train_true.shape),
        "train_keys_n": len(train_keys),
        "train_keys_singles_first": sum("+" not in key for key in train_keys[:105]),
        "train_keys_pairs_after": sum("+" in key for key in train_keys[105:]),
        "train_ridge_shape_ok": list(train_pred.shape) == [184, 5045] and train_pred.shape == train_true.shape,
        "train_key_order_ok": (
            len(train_keys) == 184
            and all("+" not in key for key in train_keys[:105])
            and all("+" in key for key in train_keys[105:])
        ),
        "chosen_alpha": float(np.asarray(predictions["chosen_alpha"]).reshape(-1)[0]),
        "test_consumed_for_final_evaluation": True,
        "do_not_reuse_26_test_pairs_for_tuning": True,
        "next_model_requires_nested_or_external_holdout": True,
        "splits": {},
    }
    checks["train_ridge_macro_mse_recomputed"] = macro_mse(train_pred, train_true)
    checks["train_ridge_macro_mse_saved"] = saved_metrics["train"]["ridge_multihot_main_effects"]["macro_mse"]
    checks["train_ridge_macro_mse_matches"] = _close(
        checks["train_ridge_macro_mse_saved"],
        checks["train_ridge_macro_mse_recomputed"],
    )

    split_map = {
        "validation": {
            "truth": "val_true",
            "control_only": "val_pred_control",
            "train_condition_balanced_perturbed_mean": "val_pred_perturbed_mean",
            "matching_single_mean": "val_pred_matching",
            "additive": "val_pred_additive",
            "ridge_multihot_main_effects": "val_pred_ridge",
        },
        "test": {
            "truth": "test_true",
            "control_only": "test_pred_control",
            "train_condition_balanced_perturbed_mean": "test_pred_perturbed_mean",
            "matching_single_mean": "test_pred_matching",
            "additive": "test_pred_additive",
            "ridge_multihot_main_effects": "test_pred_ridge",
        },
    }
    for split, names in split_map.items():
        truth = np.asarray(predictions[names["truth"]], dtype=np.float64)
        split_out = {}
        for model, key in names.items():
            if model == "truth":
                continue
            recomputed = recompute_split_metrics(
                np.asarray(predictions[key], dtype=np.float64),
                truth,
            )
            saved = saved_metrics[split][model]
            split_out[model] = {
                "saved_macro_mse": saved["macro_mse"],
                "recomputed_macro_mse": recomputed["macro_mse"],
                "macro_mse_matches": _close(saved["macro_mse"], recomputed["macro_mse"]),
                "saved_macro_mae": saved["macro_mae"],
                "recomputed_macro_mae": recomputed["macro_mae"],
                "macro_mae_matches": _close(saved["macro_mae"], recomputed["macro_mae"]),
                "saved_macro_delta_pearson": saved.get("macro_delta_pearson"),
                "recomputed_macro_delta_pearson": recomputed["macro_delta_pearson"],
                "macro_delta_pearson_matches": _close(
                    saved.get("macro_delta_pearson"),
                    recomputed["macro_delta_pearson"],
                ),
                "n_conditions": recomputed["n_conditions"],
            }
        checks["splits"][split] = split_out
    checks["all_checked_metrics_match"] = bool(
        checks["train_ridge_macro_mse_matches"]
        and checks["train_ridge_shape_ok"]
        and checks["train_key_order_ok"]
        and all(
            item["macro_mse_matches"] and item["macro_mae_matches"]
            for split in checks["splits"].values()
            for item in split.values()
        )
    )
    return checks


def write_diagnostics_report(
    own_counts: dict,
    own_frame: pd.DataFrame,
    shift_summary: dict,
    projection: pd.DataFrame,
    orientation_singles: pd.DataFrame,
    orientation_pairs: pd.DataFrame,
    rankings: pd.DataFrame,
    example_meta: list[dict],
    examples: pd.DataFrame,
    neighbors: pd.DataFrame,
    svd_signed: pd.DataFrame,
    metric_checks: dict,
    go_skip: dict,
    pair_summary: pd.DataFrame,
) -> None:
    pos_examples = own_frame.loc[own_frame["sign"] == "positive"].sort_values("effect", ascending=False).head(3)
    neg_examples = own_frame.loc[own_frame["sign"] == "negative"].sort_values("effect").head(3)
    top_proj = projection.head(5)
    largest = rankings.loc[rankings["list"] == "largest_residual_rms"]
    lowest = rankings.loc[rankings["list"] == "lowest_residual_over_response_rms"]
    if len(orientation_pairs) and "label_a" in orientation_pairs.columns:
        range_rows = orientation_pairs.loc[orientation_pairs["label_a"] == "RANGE"]
    else:
        range_rows = orientation_pairs.iloc[0:0]

    own_table = _md_table(
        [
            {
                "target": row.target_gene,
                "effect": row.effect,
                "se": row.conditional_cell_sampling_se,
                "lane+": row.positive_lane_fraction,
                "sign": row.sign,
            }
            for row in pd.concat([pos_examples, neg_examples]).itertuples(index=False)
        ],
        ["target", "effect", "se", "lane+", "sign"],
    )
    proj_table = _md_table(
        [
            {
                "pair": row.canonical_perturbation,
                "frac": row.squared_norm_fraction,
                "cos_to_-c": row.cosine_to_neg_shared_shift,
                "ortho_rms": row.orthogonal_residual_rms,
                "resid_rms": row.residual_rms,
            }
            for row in top_proj.itertuples(index=False)
        ],
        ["pair", "frac", "cos_to_-c", "ortho_rms", "resid_rms"],
    )
    large_table = _md_table(
        [
            {
                "rank": int(row.rank),
                "pair": row.canonical_perturbation,
                "residual_rms": row.residual_rms,
                "observed_rms": row.observed_rms,
                "resid/response": row.residual_over_response_rms,
            }
            for row in largest.itertuples(index=False)
        ],
        ["rank", "pair", "residual_rms", "observed_rms", "resid/response"],
    )
    low_table = _md_table(
        [
            {
                "rank": int(row.rank),
                "pair": row.canonical_perturbation,
                "residual_rms": row.residual_rms,
                "observed_rms": row.observed_rms,
                "resid/response": row.residual_over_response_rms,
            }
            for row in lowest.itertuples(index=False)
        ],
        ["rank", "pair", "residual_rms", "observed_rms", "resid/response"],
    )
    neighbor_table = _md_table(
        [
            {
                "rank": int(row.rank),
                "gene_a": row.gene_a,
                "gene_b": row.gene_b,
                "cosine": row.cosine,
            }
            for row in neighbors.itertuples(index=False)
        ],
        ["rank", "gene_a", "gene_b", "cosine"],
    )
    ori_table = _md_table(
        [
            {
                "gene": row.canonical_perturbation,
                "profile_rms": row.profile_rms,
                "pooled_rms": row.pooled_single_effect_rms,
                "ratio": row.profile_over_pooled_single_rms,
            }
            for row in orientation_singles.sort_values(
                "profile_over_pooled_single_rms",
                ascending=False,
            ).head(8).itertuples(index=False)
        ],
        ["gene", "profile_rms", "pooled_rms", "ratio"],
    ) if len(orientation_singles) else "No dual-orientation singles were available."
    range_table = _md_table(
        [
            {
                "pair": row.canonical_perturbation,
                "min": row.residual_rms_min,
                "max": row.residual_rms_max,
                "canonical": row.canonical_residual_rms,
                "n": int(row.n_combinations),
            }
            for row in range_rows.sort_values("residual_rms_max", ascending=False).head(8).itertuples(index=False)
        ],
        ["pair", "min", "max", "canonical", "n"],
    ) if len(range_rows) else "No pair had an orientation alternative."

    example_blocks = []
    for slot in example_meta:
        key = slot["canonical_perturbation"]
        block = examples.loc[examples["canonical_perturbation"] == key]
        table = _md_table(
            [
                {
                    "gene": row.gene_name,
                    "A": row.effect_a,
                    "B": row.effect_b,
                    "AB": row.observed_ab,
                    "A+B": row.additive,
                    "residual": row.residual,
                    "CI": f"[{_fmt(row.ci95_low)}, {_fmt(row.ci95_high)}]",
                    "sign_frac": row.same_sign_lane_fraction,
                }
                for row in block.itertuples(index=False)
            ],
            ["gene", "A", "B", "AB", "A+B", "residual", "CI", "sign_frac"],
        )
        example_blocks.append(
            f"### {slot['role']}: `{key}`\n\n"
            f"Residual RMS {_fmt(slot['residual_rms'])}; observed RMS {_fmt(slot['observed_rms'])}; "
            f"residual/response {_fmt(slot['residual_over_response_rms'])}. "
            "CI is the conditional cell-sampling interval on the interaction residual. No p-value.\n\n"
            f"{table}"
        )
    svd_lines = []
    for row in svd_signed.itertuples(index=False):
        svd_lines.append(
            f"- PC{int(row.component)} {row.pole}: {row.genes}"
        )

    match = "PASS" if metric_checks.get("all_checked_metrics_match") else "FAIL"
    text = f"""# Modeling diagnostics

Generated by `scripts/modeling_diagnostics.py` from frozen `reports/modeling/` artifacts.
This is a **descriptive** follow-up. It is not a new test metric, not a corrected predictor, and not a reason to retune alpha 0.1.

The 26 B-protocol test pairs are **consumed**. Next models need nested resampling or new external held-out data. Do not reuse those 26 pairs for tuning.

## 1. Own-target single effects

Of {own_counts["n_mapped"]} / {own_counts["n_targets"]} measured targets, **{own_counts["n_positive"]}** own-target effects are positive and **{own_counts["n_negative"]}** are negative ({own_counts["n_zero"]} zero; {own_counts["n_missing"]} unmapped). These are normalized log1p own-locus changes with conditional cell-sampling SE/CI. They are **not** molecular validation of guide activity, protein overexpression, or signature specificity.

{own_table}

Full table: [modeling/diagnostics/own_target_single_effects.csv](modeling/diagnostics/own_target_single_effects.csv).

## 2. Shared-shift alignment of pair residuals

`c` is the unweighted mean of the 105 training single-gene delta centroids versus **train** `mu0` from `predictions.npz`. Pair residuals `r` are the all-data lane-balanced interactions versus all-data control. Projection onto `c` is a descriptive alignment check. A high squared-norm fraction is **not** confounding proof.

Median squared-norm fraction: {_fmt(shift_summary["median_squared_norm_fraction"])}. Median cosine(`r`, `-c`): {_fmt(shift_summary["median_cosine_to_neg_shared_shift"])}. Largest fraction: `{shift_summary["largest_squared_norm_fraction_pair"]}` = {_fmt(shift_summary["largest_squared_norm_fraction"])}.

{proj_table}

Full table: [modeling/diagnostics/shared_shift_projections.csv](modeling/diagnostics/shared_shift_projections.csv).

## 3. Orientation and pair rankings

{len(orientation_singles)} dual-orientation singles: RMS between the two original within-lane mean profiles, as a ratio to the pooled single-effect RMS. Descriptive only.

{ori_table}

Alternative additive predictions use each original single label (at most four combinations). Residual RMS range is reported next to the canonical residual. No favorable orientation is selected.

{range_table}

Largest residual RMS (do not ignore observed RMS):

{large_table}

Lowest residual/response ratios, with observed RMS so weak responses are not called additive:

{low_table}

## 4. Biological examples

At most three descriptive pairs. Uncertainty is not a p-value.

{chr(10).join(example_blocks) if example_blocks else "No example pairs were selected."}

## 5. Neighboring single signatures and SVD modules

Top 10 off-diagonal cosine pairs among the 105 single-effect signatures:

{neighbor_table}

Existing SVD top-gene modules, grouped by loading sign. Data-only labels; no named pathway enrichment.

{chr(10).join(svd_lines)}

Local GO graph ({_fmt(go_skip.get("edges"))} edges, {_fmt(go_skip.get("nodes"))} nodes; columns {go_skip.get("columns")}) has no named gene-set terms. Named GO enrichment is **skipped**.

## 6. Independent metric verification

Saved val/test scorecards were recomputed from `predictions.npz` with array operations only. Ridge was **not** refit. Train ridge shape/order check: shape {metric_checks["train_ridge_shape"]}, singles-first {metric_checks["train_keys_singles_first"]}, pairs-after {metric_checks["train_keys_pairs_after"]}. Overall match: **{match}**.

Details: [modeling/diagnostics/metric_verification.json](modeling/diagnostics/metric_verification.json).

## Unknown

Whether residuals that align with `c` are construct, composition, or biology; whether orientation RMS ratios are experimental nonequivalence; whether own-target sign is on-target CRISPRa. Those questions are not answered by this diagnostic.
"""
    REPORT_PATH.write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    print("diagnostics: reading frozen modeling artifacts", flush=True)
    run_diagnostics()
    print("diagnostics: wrote " + str(DIAG_DIR), flush=True)
    print("diagnostics: wrote " + str(REPORT_PATH), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
