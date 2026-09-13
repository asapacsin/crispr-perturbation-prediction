"""CLI for authorized Package-B baselines and post-fit descriptive effects.

Stages: prepare, fit, describe, or all (fit before describe).
Loads sparse stored X only. Never loads the dense counts layer. Does not
read uns, mutate data/norman, or change historical reports or the schema.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from collections.abc import Mapping, Sequence
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

try:
    import yaml
except ImportError as exc:  # pragma: no cover - environment check
    raise SystemExit(
        "ERROR: PyYAML is required to read config/modeling_protocol.yaml. "
        "This script will not pip-install it. "
        f"Original import error: {exc}"
    ) from exc

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

from dataset_schema_lib import (
    canonicalize_condition,
    combination_split_errors,
    n_targets,
)
from modeling_lib import (
    BOOTSTRAP_SEED,
    LANE_IDS,
    N_BOOTSTRAP,
    RIDGE_ALPHAS,
    SIGNIFICANCE_NOT_ESTIMABLE,
    Z_NORMAL_95,
    additive_delta,
    condition_signature_svd,
    compute_group_lane_stats,
    cosine_similarity,
    delta_pearson,
    dual_orientation_single_keys,
    group_orientation_labels,
    json_safe,
    lane_balanced_mean,
    leave_one_lane_out_rms,
    macro_mae,
    macro_mse,
    matching_single_delta,
    multi_hot_rows,
    pair_interaction_lane,
    pair_partners,
    paired_bootstrap_mean_difference,
    per_condition_mae,
    per_condition_mse,
    run_ridge_selection,
    single_effect_lane,
    skill_score,
    summarize_lane_contrasts,
    vector_rms,
    zero_safe_ratio,
)
from validate_dataset_schema import load_obs_frame, load_var_frame, matrix_shape

DEFAULT_PROTOCOL = ROOT / "config" / "modeling_protocol.yaml"
MODELING_DIR = ROOT / "reports" / "modeling"
PREPARE_DIR = MODELING_DIR / "prepare"
FIT_DIR = MODELING_DIR / "fit"
DESCRIBE_DIR = MODELING_DIR / "describe"
FIGURE_DIR = DESCRIBE_DIR / "figures"
REPORT_PATH = ROOT / "reports" / "MODELING_REPORT.md"
LANE_STATS_PATH = PREPARE_DIR / "lane_stats.h5"

BASELINE_ORDER = (
    "control_only",
    "train_condition_balanced_perturbed_mean",
    "matching_single_mean",
    "additive",
    "ridge_multihot_main_effects",
)


class ModelingError(RuntimeError):
    pass


def load_protocol(path: Path = DEFAULT_PROTOCOL) -> dict:
    protocol = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not protocol.get("authorized"):
        raise ModelingError("modeling_protocol.yaml is not authorized")
    return protocol


def load_dataset_schema() -> dict:
    return yaml.safe_load((ROOT / "config" / "dataset_schema.yaml").read_text(encoding="utf-8"))


def load_csr_X(h5: h5py.File) -> sp.csr_matrix:
    if "X" not in h5:
        raise ModelingError("HDF5 is missing X")
    node = h5["X"]
    if isinstance(node, h5py.Group):
        if not {"data", "indices", "indptr"} <= set(node.keys()):
            raise ModelingError(
                f"X group is not a CSR sparse matrix; keys={list(node.keys())}"
            )
        data = node["data"][()]
        indices = node["indices"][()]
        indptr = node["indptr"][()]
        return sp.csr_matrix((data, indices, indptr), shape=matrix_shape(node))
    if isinstance(node, h5py.Dataset):
        return sp.csr_matrix(node[()])
    raise ModelingError(f"cannot load X of type {type(node).__name__}")


def barcode_lane(barcode: str) -> int:
    return int(str(barcode).rsplit("-", 1)[-1])


def write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(json_safe(payload), indent=2, sort_keys=True),
        encoding="utf-8",
    )


def write_lane_stats_group(container: h5py.Group, stats: dict) -> None:
    for key, payload in stats.items():
        group = container.create_group(key)
        group.create_dataset("n", data=np.asarray(payload["n"]))
        group.create_dataset("mean", data=np.asarray(payload["mean"]))
        group.create_dataset("sample_var", data=np.asarray(payload["sample_var"]))


def read_lane_stats_group(container: h5py.Group) -> dict:
    out = {}
    for key in container.keys():
        out[key] = {
            "n": container[key]["n"][()],
            "mean": container[key]["mean"][()],
            "sample_var": container[key]["sample_var"][()],
        }
    return out


def condition_mean(entry: dict) -> np.ndarray:
    return lane_balanced_mean(entry["mean"])


def split_target_sets(obs: pd.DataFrame, column: str) -> dict[str, list[str]]:
    noncontrol = obs.loc[obs["n_targets"] > 0]
    assigned = noncontrol.groupby("canonical_perturbation", sort=True)[column].nunique()
    leaked = assigned[assigned > 1]
    if len(leaked):
        raise ModelingError(f"canonical sets in multiple splits: {list(leaked.index)}")
    mapping = noncontrol.groupby("canonical_perturbation", sort=True)[column].first()
    out = {"train": [], "validation": [], "test": []}
    for key, split in mapping.items():
        if split not in out:
            raise ModelingError(f"unexpected split {split!r} for {key}")
        out[split].append(str(key))
    for split in out:
        singles = [key for key in out[split] if "+" not in key]
        pairs = [key for key in out[split] if "+" in key]
        out[split] = sorted(singles) + sorted(pairs)
    return out


def baseline_predictions(
    keys: Sequence[str],
    true_means: np.ndarray,
    mu0: np.ndarray,
    single_means: dict[str, np.ndarray],
    train_perturbed_mean: np.ndarray,
    ridge_pred: np.ndarray,
) -> dict[str, np.ndarray]:
    n_cond, n_genes = true_means.shape
    zeros = np.zeros((n_cond, n_genes), dtype=np.float64)
    common = np.broadcast_to(train_perturbed_mean - mu0, (n_cond, n_genes)).copy()
    matching = np.empty_like(true_means)
    additive = np.empty_like(true_means)
    for i, key in enumerate(keys):
        gene_a, gene_b = pair_partners(key)
        delta_a = single_means[gene_a] - mu0
        delta_b = single_means[gene_b] - mu0
        matching[i] = matching_single_delta(delta_a, delta_b)
        additive[i] = additive_delta(delta_a, delta_b)
    return {
        "control_only": zeros,
        "train_condition_balanced_perturbed_mean": common,
        "matching_single_mean": matching,
        "additive": additive,
        "ridge_multihot_main_effects": ridge_pred,
        "truth": true_means - mu0,
    }


def scorecard(pred: np.ndarray, truth: np.ndarray, keys: Sequence[str], model: str, split: str) -> pd.DataFrame:
    mse = per_condition_mse(pred, truth)
    mae = per_condition_mae(pred, truth)
    rows = []
    for i, key in enumerate(keys):
        rows.append(
            {
                "split": split,
                "model": model,
                "canonical_perturbation": key,
                "mse": float(mse[i]),
                "mae": float(mae[i]),
                "delta_pearson": delta_pearson(pred[i], truth[i]),
                "cosine": cosine_similarity(pred[i], truth[i]),
            }
        )
    return pd.DataFrame(rows)


def summary_from_scorecard(frame: pd.DataFrame, split: str) -> dict:
    sub = frame.loc[frame["split"] == split]
    out = {}
    for model, block in sub.groupby("model", sort=False):
        pearson = block["delta_pearson"].to_numpy(dtype=np.float64)
        cosine = block["cosine"].to_numpy(dtype=np.float64)
        out[model] = {
            "macro_mse": float(block["mse"].mean()),
            "macro_mae": float(block["mae"].mean()),
            "macro_delta_pearson": float(np.nanmean(pearson)) if np.any(np.isfinite(pearson)) else None,
            "macro_cosine": float(np.nanmean(cosine)) if np.any(np.isfinite(cosine)) else None,
            "n_pearson_defined": int(np.isfinite(pearson).sum()),
            "n_cosine_defined": int(np.isfinite(cosine).sum()),
            "n_conditions": int(len(block)),
        }
    models = out
    if "ridge_multihot_main_effects" in models:
        ridge_mse = models["ridge_multihot_main_effects"]["macro_mse"]
        for baseline in (
            "train_condition_balanced_perturbed_mean",
            "additive",
        ):
            if baseline in models:
                models["ridge_multihot_main_effects"][f"skill_vs_{baseline}"] = skill_score(
                    ridge_mse,
                    models[baseline]["macro_mse"],
                )
        if "additive" in models:
            models["additive"]["skill_vs_train_condition_balanced_perturbed_mean"] = skill_score(
                models["additive"]["macro_mse"],
                models["train_condition_balanced_perturbed_mean"]["macro_mse"],
            )
    return models


def stage_prepare(protocol: dict) -> dict:
    schema = load_dataset_schema()
    h5ad_path = ROOT / schema["identity"]["local_h5ad"]
    if not h5ad_path.is_file():
        raise ModelingError(f"H5AD not found: {h5ad_path}")
    PREPARE_DIR.mkdir(parents=True, exist_ok=True)
    (MODELING_DIR / "protocol_snapshot.yaml").write_text(yaml.safe_dump(protocol, sort_keys=False), encoding="utf-8")

    split_path = ROOT / protocol["split_manifest"]
    splits = pd.read_csv(split_path)
    if "cell_barcode" not in splits.columns:
        raise ModelingError("split manifest is missing cell_barcode")

    with h5py.File(h5ad_path, "r") as h5:
        if "layers" in h5 and "counts" in h5["layers"]:
            # Shape check only; do not read the counts payload.
            counts_shape = matrix_shape(h5["layers"]["counts"])
        else:
            counts_shape = None
        x_shape = matrix_shape(h5["X"])
        obs = load_obs_frame(h5)
        var = load_var_frame(h5)
        X = load_csr_X(h5)

    if x_shape != tuple(schema["identity"]["shape"]):
        raise ModelingError(f"X shape {x_shape} != schema {schema['identity']['shape']}")
    if counts_shape not in {None, tuple(schema["identity"]["shape"])}:
        raise ModelingError(f"counts shape {counts_shape} != schema")
    if X.shape != x_shape:
        raise ModelingError(f"loaded CSR shape {X.shape} != {x_shape}")
    if list(var["gene_id"]) != list(dict.fromkeys(var["gene_id"])):
        raise ModelingError("gene_id order is not unique")
    if len(var) != protocol["expression"]["n_features"]:
        raise ModelingError("feature count does not match the predeclared 5045")

    obs["canonical_perturbation"] = obs["condition"].map(canonicalize_condition)
    obs["n_targets"] = obs["condition"].map(n_targets)
    obs["barcode_group"] = obs["cell_barcode"].map(barcode_lane)
    obs = obs.merge(splits[["cell_barcode", protocol["split"]["column"]]], on="cell_barcode", how="left")
    split_col = protocol["split"]["column"]
    if obs[split_col].isna().any():
        raise ModelingError("split manifest does not cover every barcode")

    b_rows = list(zip(obs["condition"], obs[split_col].astype(str)))
    leak = combination_split_errors(b_rows)
    if leak:
        raise ModelingError("B split leakage: " + "; ".join(leak[:8]))

    sets = split_target_sets(obs, split_col)
    expected = protocol["split"]["expected"]
    for split_name, spec in expected.items():
        actual_cells = int((obs[split_col] == split_name).sum())
        if actual_cells != spec["cells"]:
            raise ModelingError(f"{split_name} has {actual_cells} cells, protocol says {spec['cells']}")
        singles = [key for key in sets[split_name] if "+" not in key]
        pairs = [key for key in sets[split_name] if "+" in key]
        if len(singles) != spec["single_target_sets"] or len(pairs) != spec["double_target_sets"]:
            raise ModelingError(
                f"{split_name} sets singles={len(singles)} pairs={len(pairs)} "
                f"!= protocol {spec['single_target_sets']}/{spec['double_target_sets']}"
            )

    grouped = group_orientation_labels(obs["condition"].tolist())
    dual = dual_orientation_single_keys(grouped)
    if len(dual) != protocol["orientation_sensitivity"]["expected_dual_orientation_singles"]:
        raise ModelingError(
            f"dual-orientation singles {len(dual)} != "
            f"{protocol['orientation_sensitivity']['expected_dual_orientation_singles']}"
        )

    print("prepare: computing canonical lane stats", flush=True)
    canonical_stats = compute_group_lane_stats(
        X,
        obs["canonical_perturbation"].tolist(),
        obs["barcode_group"].tolist(),
    )
    print("prepare: computing original-label lane stats", flush=True)
    original_stats = compute_group_lane_stats(
        X,
        obs["condition"].astype(str).tolist(),
        obs["barcode_group"].tolist(),
    )
    split_controls = {}
    for split_name in ("train", "validation", "test"):
        mask = (obs["canonical_perturbation"] == "ctrl") & (obs[split_col] == split_name)
        split_controls[split_name] = compute_group_lane_stats(
            X[mask.to_numpy()],
            obs.loc[mask, "canonical_perturbation"].tolist(),
            obs.loc[mask, "barcode_group"].tolist(),
        )

    with h5py.File(LANE_STATS_PATH, "w") as store:
        write_lane_stats_group(store.create_group("canonical"), canonical_stats)
        write_lane_stats_group(store.create_group("original"), original_stats)
        split_group = store.create_group("split_control")
        for split_name, stats in split_controls.items():
            write_lane_stats_group(split_group.create_group(split_name), stats)

    var_out = var.copy()
    var_out["gene_index"] = np.arange(len(var_out), dtype=np.int64)
    var_out.to_csv(PREPARE_DIR / "genes.csv", index=False, encoding="utf-8")

    examples = obs.iloc[:12].loc[
        :,
        [
            "cell_barcode",
            "condition",
            "canonical_perturbation",
            "n_targets",
            "barcode_group",
            split_col,
            "control",
            "cell_type",
            "dose_val",
            "condition_name",
        ],
    ].copy()
    examples["stored_cell_type_is_dummy"] = True
    examples.to_csv(PREPARE_DIR / "cell_examples.csv", index=False, encoding="utf-8")

    size_rows = []
    for key, block in obs.groupby("canonical_perturbation", sort=True):
        lane_counts = block["barcode_group"].value_counts().to_dict()
        originals = sorted(block["condition"].astype(str).unique())
        size_rows.append(
            {
                "canonical_perturbation": key,
                "n_targets": int(n_targets(key)) if key != "ctrl" else 0,
                "original_condition_labels": "|".join(originals),
                "n_original_labels": len(originals),
                "total_n": int(len(block)),
                "B_heldout_combination": (
                    "control_split_by_cell"
                    if key == "ctrl"
                    else block[split_col].iloc[0]
                ),
                **{f"n_lane_{lane}": int(lane_counts.get(lane, 0)) for lane in LANE_IDS},
            }
        )
    pd.DataFrame(size_rows).to_csv(
        PREPARE_DIR / "condition_sizes.csv",
        index=False,
        encoding="utf-8",
    )
    pd.DataFrame(
        {
            "split": [split for split, keys in sets.items() for _ in keys],
            "canonical_perturbation": [key for keys in sets.values() for key in keys],
            "n_targets": [
                1 if "+" not in key else 2
                for keys in sets.values()
                for key in keys
            ],
        }
    ).to_csv(PREPARE_DIR / "split_conditions.csv", index=False, encoding="utf-8")

    meta = {
        "h5ad": str(h5ad_path.relative_to(ROOT)).replace("\\", "/"),
        "n_obs": int(len(obs)),
        "n_vars": int(len(var)),
        "x_shape": list(X.shape),
        "x_nnz": int(X.nnz),
        "counts_payload_loaded": False,
        "uns_read": False,
        "split_column": split_col,
        "split_seed": protocol["split"]["seed"],
        "n_original_condition_labels": int(obs["condition"].nunique()),
        "n_canonical_target_sets_including_control": int(obs["canonical_perturbation"].nunique()),
        "n_control_cells": int((obs["canonical_perturbation"] == "ctrl").sum()),
        "n_single_cells": int((obs["n_targets"] == 1).sum()),
        "n_double_cells": int((obs["n_targets"] == 2).sum()),
        "dual_orientation_singles": dual,
        "split_sets": sets,
        "gene_id_head": var["gene_id"].head(5).tolist(),
        "loaded_via": "h5py_csr_X_only",
    }
    write_json(PREPARE_DIR / "prepare_meta.json", meta)
    print("prepare: wrote " + str(PREPARE_DIR), flush=True)
    return meta


def _require_prepare() -> None:
    if not LANE_STATS_PATH.is_file() or not (PREPARE_DIR / "prepare_meta.json").is_file():
        raise ModelingError("run prepare before fit")


def _require_fit() -> None:
    if not (FIT_DIR / "chosen_alpha.json").is_file():
        raise ModelingError("run fit before describe; describe cannot select alpha")


def stage_fit(protocol: dict) -> dict:
    _require_prepare()
    FIT_DIR.mkdir(parents=True, exist_ok=True)
    schema = load_dataset_schema()
    target_names = list(schema["targets"]["names"])
    if len(target_names) != 105:
        raise ModelingError("schema target list is not 105 names")

    genes = pd.read_csv(PREPARE_DIR / "genes.csv")
    gene_ids = genes["gene_id"].tolist()
    prepare_meta = json.loads((PREPARE_DIR / "prepare_meta.json").read_text(encoding="utf-8"))
    sets = prepare_meta["split_sets"]
    train_keys = list(sets["train"])
    val_keys = list(sets["validation"])
    test_keys = list(sets["test"])
    if len(train_keys) != 184 or len(val_keys) != 26 or len(test_keys) != 26:
        raise ModelingError("prepared split set counts do not match 184/26/26")
    if any("+" not in key for key in val_keys + test_keys):
        raise ModelingError("validation/test must be pairs only")

    with h5py.File(LANE_STATS_PATH, "r") as store:
        canonical = read_lane_stats_group(store["canonical"])
        train_ctrl = read_lane_stats_group(store["split_control"]["train"])
        val_ctrl = read_lane_stats_group(store["split_control"]["validation"])
        test_ctrl = read_lane_stats_group(store["split_control"]["test"])

    mu0 = condition_mean(train_ctrl["ctrl"])
    val_ctrl_mean = condition_mean(val_ctrl["ctrl"])
    test_ctrl_mean = condition_mean(test_ctrl["ctrl"])
    single_means = {key: condition_mean(canonical[key]) for key in target_names}
    train_means = np.stack([condition_mean(canonical[key]) for key in train_keys])
    val_means = np.stack([condition_mean(canonical[key]) for key in val_keys])
    test_means = np.stack([condition_mean(canonical[key]) for key in test_keys])
    train_deltas = train_means - mu0
    train_perturbed_mean = train_means.mean(axis=0)

    train_x = multi_hot_rows(train_keys, target_names)
    val_x = multi_hot_rows(val_keys, target_names)
    test_x = multi_hot_rows(test_keys, target_names)
    if np.any(train_x.sum(axis=0) == 0):
        raise ModelingError("a target gene is absent from the training multi-hot design")
    if train_x.shape != (184, 105):
        raise ModelingError(f"train design {train_x.shape} != (184, 105)")

    ridge = run_ridge_selection(
        train_x,
        train_deltas,
        val_x,
        val_means - mu0,
        alphas=tuple(protocol["ridge"]["alphas"]),
        test_x=test_x,
        test_y=np.full(test_means.shape, np.nan),
    )
    if ridge["alpha"] not in set(protocol["ridge"]["alphas"]) and ridge["alpha"] not in set(
        float(a) for a in protocol["ridge"]["alphas"]
    ):
        raise ModelingError("selected alpha is not on the predeclared grid")

    val_pred = baseline_predictions(
        val_keys,
        val_means,
        mu0,
        single_means,
        train_perturbed_mean,
        ridge["validation_prediction"],
    )
    test_pred = baseline_predictions(
        test_keys,
        test_means,
        mu0,
        single_means,
        train_perturbed_mean,
        ridge["test_prediction"],
    )
    train_ridge = ridge["train_prediction"]
    train_pred = {
        "control_only": np.zeros_like(train_deltas),
        "train_condition_balanced_perturbed_mean": np.broadcast_to(
            train_perturbed_mean - mu0,
            train_deltas.shape,
        ).copy(),
        "matching_single_mean": np.stack(
            [
                matching_single_delta(single_means[pair_partners(key)[0]] - mu0, single_means[pair_partners(key)[1]] - mu0)
                if "+" in key
                else single_means[key] - mu0
                for key in train_keys
            ]
        ),
        "additive": np.stack(
            [
                additive_delta(single_means[pair_partners(key)[0]] - mu0, single_means[pair_partners(key)[1]] - mu0)
                if "+" in key
                else single_means[key] - mu0
                for key in train_keys
            ]
        ),
        "ridge_multihot_main_effects": train_ridge,
        "truth": train_deltas,
    }

    metric_frames = []
    for split, payload, keys in (
        ("train", train_pred, train_keys),
        ("validation", val_pred, val_keys),
        ("test", test_pred, test_keys),
    ):
        truth = payload["truth"]
        for model in BASELINE_ORDER:
            metric_frames.append(scorecard(payload[model], truth, keys, model, split))
    metrics = pd.concat(metric_frames, ignore_index=True)
    metrics.to_csv(FIT_DIR / "per_condition_metrics.csv", index=False, encoding="utf-8")

    summaries = {
        split: summary_from_scorecard(metrics, split)
        for split in ("train", "validation", "test")
    }
    test_ridge = metrics[
        (metrics["split"] == "test")
        & (metrics["model"] == "ridge_multihot_main_effects")
    ].sort_values("canonical_perturbation")
    test_add = metrics[
        (metrics["split"] == "test") & (metrics["model"] == "additive")
    ].sort_values("canonical_perturbation")
    if list(test_ridge["canonical_perturbation"]) != list(test_add["canonical_perturbation"]):
        raise ModelingError("test ridge/additive condition order mismatch")
    bootstrap = paired_bootstrap_mean_difference(
        test_ridge["mse"].to_numpy(),
        test_add["mse"].to_numpy(),
        seed=protocol["metrics"]["bootstrap"]["seed"],
        n_samples=protocol["metrics"]["bootstrap"]["n_samples"],
    )

    pd.DataFrame(ridge["validation_curve"]).to_csv(
        FIT_DIR / "validation_curve.csv",
        index=False,
        encoding="utf-8",
    )
    for name, keys in (
        ("train_conditions.csv", train_keys),
        ("val_conditions.csv", val_keys),
        ("test_conditions.csv", test_keys),
    ):
        pd.DataFrame({"canonical_perturbation": keys}).to_csv(
            FIT_DIR / name,
            index=False,
            encoding="utf-8",
        )

    np.savez_compressed(
        FIT_DIR / "predictions.npz",
        gene_ids=np.asarray(gene_ids),
        target_names=np.asarray(target_names),
        mu0=mu0,
        train_keys=np.asarray(train_keys),
        val_keys=np.asarray(val_keys),
        test_keys=np.asarray(test_keys),
        train_true=train_deltas,
        val_true=val_pred["truth"],
        test_true=test_pred["truth"],
        train_pred_ridge=train_ridge,
        val_pred_ridge=val_pred["ridge_multihot_main_effects"],
        test_pred_ridge=test_pred["ridge_multihot_main_effects"],
        val_pred_additive=val_pred["additive"],
        test_pred_additive=test_pred["additive"],
        val_pred_matching=val_pred["matching_single_mean"],
        test_pred_matching=test_pred["matching_single_mean"],
        val_pred_perturbed_mean=val_pred["train_condition_balanced_perturbed_mean"],
        test_pred_perturbed_mean=test_pred["train_condition_balanced_perturbed_mean"],
        val_pred_control=val_pred["control_only"],
        test_pred_control=test_pred["control_only"],
        ridge_weights=ridge["weights"],
        ridge_intercept=ridge["intercept"],
        chosen_alpha=np.array([ridge["alpha"]]),
        train_perturbed_mean=train_perturbed_mean,
    )

    chosen = {
        "alpha": ridge["alpha"],
        "selection_metric": "validation_macro_mse",
        "validation_curve": ridge["validation_curve"],
        "refit_train_plus_val": False,
        "test_evaluated_once": True,
        "test_used_in_selection": False,
        "describe_used_in_fit": False,
        "n_train_conditions": 184,
        "n_features": 105,
        "encoding": "multi_hot_main_effects_only",
        "intercept": "unpenalized",
        "center_design_and_response": True,
        "val_test_controls_used_as_reference": False,
        "val_control_vs_train_mu0": {
            "cosine": cosine_similarity(val_ctrl_mean, mu0),
            "delta_pearson": delta_pearson(val_ctrl_mean, mu0),
            "mse": float(np.mean((val_ctrl_mean - mu0) ** 2)),
            "role": "diagnostic_only",
        },
        "test_control_vs_train_mu0": {
            "cosine": cosine_similarity(test_ctrl_mean, mu0),
            "delta_pearson": delta_pearson(test_ctrl_mean, mu0),
            "mse": float(np.mean((test_ctrl_mean - mu0) ** 2)),
            "role": "diagnostic_only",
        },
        "same_reference_mse_identity": (
            "MSE(mu-r, hatmu-r) equals MSE(mu, hatmu) for shared r; "
            "reported numbers are on deltas versus train mu0"
        ),
    }
    write_json(FIT_DIR / "chosen_alpha.json", chosen)
    write_json(FIT_DIR / "summary_metrics.json", summaries)
    write_json(FIT_DIR / "bootstrap_vs_additive.json", bootstrap)
    write_json(
        FIT_DIR / "fit_meta.json",
        {
            "split_column": protocol["split"]["column"],
            "split_seed": protocol["split"]["seed"],
            "ridge_alphas": list(protocol["ridge"]["alphas"]),
            "chosen_alpha": ridge["alpha"],
            "bootstrap_seed": BOOTSTRAP_SEED,
            "n_bootstrap": N_BOOTSTRAP,
            "primary_metric": "macro_condition_mse",
            "all_data_describe_used": False,
        },
    )
    print(f"fit: selected alpha={ridge['alpha']}", flush=True)
    return chosen


def _long_frame(
    keys: Sequence[str],
    genes: pd.DataFrame,
    originals: Mapping[str, Sequence[str]],
    means: np.ndarray,
    effects: np.ndarray,
    se: np.ndarray,
    ci_low: np.ndarray,
    ci_high: np.ndarray,
    lane_frac: np.ndarray,
    lane_frac_name: str,
    min_n: Sequence[int],
    total_n: Sequence[int],
    n_lanes_used: Sequence[int],
    se_flag: Sequence[bool],
    extra_cols: Sequence[tuple[str, np.ndarray]] | None = None,
) -> pd.DataFrame:
    n_p, n_g = effects.shape
    frame = {
        "canonical_perturbation": np.repeat(np.asarray(keys, dtype=object), n_g),
        "original_condition_labels": np.repeat(
            np.asarray(["|".join(originals[key]) for key in keys], dtype=object),
            n_g,
        ),
        "gene_id": np.tile(genes["gene_id"].to_numpy(), n_p),
        "gene_name": np.tile(genes["gene_name"].to_numpy(), n_p),
        "gene_index": np.tile(genes["gene_index"].to_numpy(), n_p),
        "condition_mean": means.ravel(),
        "effect": effects.ravel(),
        "conditional_cell_sampling_se": se.ravel(),
        "ci95_low": ci_low.ravel(),
        "ci95_high": ci_high.ravel(),
        lane_frac_name: lane_frac.ravel(),
        "minimum_lane_n": np.repeat(np.asarray(min_n, dtype=np.int64), n_g),
        "total_n": np.repeat(np.asarray(total_n, dtype=np.int64), n_g),
        "n_lanes_used": np.repeat(np.asarray(n_lanes_used, dtype=np.int64), n_g),
        "se_undefined_lane_n_lt_2": np.repeat(np.asarray(se_flag, dtype=bool), n_g),
        "significance": SIGNIFICANCE_NOT_ESTIMABLE,
    }
    if extra_cols:
        for name, values in extra_cols:
            frame[name] = values.ravel()
    return pd.DataFrame(frame)


def _try_pyplot():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        return plt
    except ImportError:
        return None


PAIR_RANKING_CAPTION = (
    "Distribution of log-space additive residual RMS across 131 pairs, "
    "and the 20 largest residuals."
)
PAIR_RANKING_UNITS = "RMS over 5,045 genes of the lane-averaged interaction residual."
PAIR_RANKING_CAVEAT = (
    "Descriptive all-data ranking. Uncertainty is conditional cell-sampling, "
    "not independent biological replication."
)


def _save_figure(fig, stem: str, caption: str, units: str, caveat: str) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    png = FIGURE_DIR / f"{stem}.png"
    fig.savefig(png, dpi=150, bbox_inches="tight")
    (FIGURE_DIR / f"{stem}.caption.md").write_text(
        f"**Caption.** {caption}\n\n**Units.** {units}\n\n**Caveat.** {caveat}\n",
        encoding="utf-8",
    )


def plot_pair_error_distribution_ranking(pair_summary: pd.DataFrame, plt):
    """Two-panel residual ranking. Constrained layout keeps y-labels off the histogram."""
    frame = pair_summary.sort_values("additive_rmse", ascending=False)
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(12.8, 5.6),
        constrained_layout=True,
        gridspec_kw={"width_ratios": [1.0, 1.55]},
    )
    axes[0].hist(frame["additive_rmse"], bins=20, color="#4C72B0")
    axes[0].set_xlabel("Additive residual RMS", fontsize=9)
    axes[0].set_ylabel("Number of pairs", fontsize=9)
    axes[0].tick_params(labelsize=8)
    ranked = frame.head(20)
    axes[1].barh(
        ranked["canonical_perturbation"][::-1],
        ranked["additive_rmse"][::-1],
        color="#C44E52",
    )
    axes[1].set_xlabel("Additive residual RMS", fontsize=9)
    axes[1].tick_params(axis="y", labelsize=6)
    axes[1].tick_params(axis="x", labelsize=8)
    fig.suptitle("Pair non-additivity ranking (all 131 pairs, descriptive)", fontsize=11)
    return fig


def refresh_pair_error_distribution_ranking() -> Path:
    """Regenerate only the ranking figure from saved pair_summaries.csv."""
    plt = _try_pyplot()
    if plt is None:
        raise ModelingError("matplotlib is required to refresh pair_error_distribution_ranking")
    summary_path = DESCRIBE_DIR / "pair_summaries.csv"
    if not summary_path.is_file():
        raise ModelingError(f"missing {summary_path}")
    fig = plot_pair_error_distribution_ranking(pd.read_csv(summary_path), plt)
    _save_figure(
        fig,
        "pair_error_distribution_ranking",
        PAIR_RANKING_CAPTION,
        PAIR_RANKING_UNITS,
        PAIR_RANKING_CAVEAT,
    )
    plt.close(fig)
    return FIGURE_DIR / "pair_error_distribution_ranking.png"


def stage_describe(protocol: dict) -> dict:
    _require_prepare()
    _require_fit()
    DESCRIBE_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    chosen = json.loads((FIT_DIR / "chosen_alpha.json").read_text(encoding="utf-8"))
    fit_summaries = json.loads((FIT_DIR / "summary_metrics.json").read_text(encoding="utf-8"))
    bootstrap = json.loads((FIT_DIR / "bootstrap_vs_additive.json").read_text(encoding="utf-8"))
    prepare_meta = json.loads((PREPARE_DIR / "prepare_meta.json").read_text(encoding="utf-8"))
    sizes = pd.read_csv(PREPARE_DIR / "condition_sizes.csv")
    genes = pd.read_csv(PREPARE_DIR / "genes.csv")
    examples = pd.read_csv(PREPARE_DIR / "cell_examples.csv")
    metrics = pd.read_csv(FIT_DIR / "per_condition_metrics.csv")
    curve = pd.read_csv(FIT_DIR / "validation_curve.csv")

    with h5py.File(LANE_STATS_PATH, "r") as store:
        canonical = read_lane_stats_group(store["canonical"])
        original = read_lane_stats_group(store["original"])

    ctrl = canonical["ctrl"]
    ctrl_mean = condition_mean(ctrl)
    originals = {
        row.canonical_perturbation: str(row.original_condition_labels).split("|")
        for row in sizes.itertuples(index=False)
    }
    size_n = {
        row.canonical_perturbation: int(row.total_n) for row in sizes.itertuples(index=False)
    }

    single_keys = [
        key for key in sorted(canonical) if key != "ctrl" and "+" not in key
    ]
    pair_keys = [key for key in sorted(canonical) if "+" in key]
    if len(single_keys) != 105 or len(pair_keys) != 131:
        raise ModelingError(f"expected 105/131, got {len(single_keys)}/{len(pair_keys)}")

    n_genes = len(genes)
    single_mean = np.empty((105, n_genes), dtype=np.float64)
    single_effect = np.empty((105, n_genes), dtype=np.float64)
    single_se = np.empty((105, n_genes), dtype=np.float64)
    single_lo = np.empty((105, n_genes), dtype=np.float64)
    single_hi = np.empty((105, n_genes), dtype=np.float64)
    single_pos = np.empty((105, n_genes), dtype=np.float64)
    single_min_n = []
    single_used = []
    single_flag = []
    single_lolo = []
    single_summary_rows = []

    for i, key in enumerate(single_keys):
        entry = canonical[key]
        contrasts = single_effect_lane(entry["mean"], ctrl["mean"])
        packed = summarize_lane_contrasts(
            contrasts,
            [entry["n"], ctrl["n"]],
            [entry["sample_var"], ctrl["sample_var"]],
            z=Z_NORMAL_95,
        )
        lolo = leave_one_lane_out_rms(contrasts)
        single_mean[i] = condition_mean(entry)
        single_effect[i] = packed["estimate"]
        single_se[i] = packed["se"]
        single_lo[i] = packed["ci95_low"]
        single_hi[i] = packed["ci95_high"]
        single_pos[i] = packed["positive_lane_fraction"]
        single_min_n.append(packed["minimum_lane_n"])
        single_used.append(packed["n_lanes_used"])
        single_flag.append(packed["se_undefined_lane_n_lt_2"])
        single_lolo.append(
            {
                "canonical_perturbation": key,
                "max_rms_change": lolo["max_rms_change"],
                "min_lolo_effect_rms": lolo["min_lolo_effect_rms"],
                "max_lolo_effect_rms": lolo["max_lolo_effect_rms"],
                "full_effect_rms": lolo["full_effect_rms"],
            }
        )
        order = np.argsort(-single_effect[i])
        top_up = order[:10]
        top_down = order[-10:][::-1]
        single_summary_rows.append(
            {
                "canonical_perturbation": key,
                "original_condition_labels": "|".join(originals[key]),
                "construct_pooling_caveat": len(originals[key]) > 1,
                "total_n": size_n[key],
                "minimum_lane_n": packed["minimum_lane_n"],
                "effect_rms": vector_rms(packed["estimate"]),
                "mean_positive_lane_fraction": float(np.mean(packed["positive_lane_fraction"])),
                "se_undefined_lane_n_lt_2": packed["se_undefined_lane_n_lt_2"],
                "lolo_max_rms_change": lolo["max_rms_change"],
                "top10_up_genes": "|".join(genes.loc[top_up, "gene_name"].astype(str)),
                "top10_down_genes": "|".join(genes.loc[top_down, "gene_name"].astype(str)),
            }
        )

    pair_mean = np.empty((131, n_genes), dtype=np.float64)
    pair_obs = np.empty((131, n_genes), dtype=np.float64)
    pair_add = np.empty((131, n_genes), dtype=np.float64)
    pair_resid = np.empty((131, n_genes), dtype=np.float64)
    pair_se = np.empty((131, n_genes), dtype=np.float64)
    pair_lo = np.empty((131, n_genes), dtype=np.float64)
    pair_hi = np.empty((131, n_genes), dtype=np.float64)
    pair_sign = np.empty((131, n_genes), dtype=np.float64)
    pair_min_n = []
    pair_used = []
    pair_flag = []
    pair_lolo = []
    pair_summary_rows = []
    effect_by_single = {key: single_effect[i] for i, key in enumerate(single_keys)}

    for i, key in enumerate(pair_keys):
        gene_a, gene_b = pair_partners(key)
        entry = canonical[key]
        a_entry = canonical[gene_a]
        b_entry = canonical[gene_b]
        contrasts = pair_interaction_lane(
            entry["mean"],
            a_entry["mean"],
            b_entry["mean"],
            ctrl["mean"],
        )
        packed = summarize_lane_contrasts(
            contrasts,
            [entry["n"], a_entry["n"], b_entry["n"], ctrl["n"]],
            [
                entry["sample_var"],
                a_entry["sample_var"],
                b_entry["sample_var"],
                ctrl["sample_var"],
            ],
            z=Z_NORMAL_95,
        )
        lolo = leave_one_lane_out_rms(contrasts)
        mean_vec = condition_mean(entry)
        observed = mean_vec - ctrl_mean
        additive = additive_delta(effect_by_single[gene_a], effect_by_single[gene_b])
        residual = packed["estimate"]
        pair_mean[i] = mean_vec
        pair_obs[i] = observed
        pair_add[i] = additive
        pair_resid[i] = residual
        pair_se[i] = packed["se"]
        pair_lo[i] = packed["ci95_low"]
        pair_hi[i] = packed["ci95_high"]
        pair_sign[i] = packed["same_sign_lane_fraction"]
        pair_min_n.append(packed["minimum_lane_n"])
        pair_used.append(packed["n_lanes_used"])
        pair_flag.append(packed["se_undefined_lane_n_lt_2"])
        pair_lolo.append(
            {
                "canonical_perturbation": key,
                "max_rms_change": lolo["max_rms_change"],
                "min_lolo_effect_rms": lolo["min_lolo_effect_rms"],
                "max_lolo_effect_rms": lolo["max_lolo_effect_rms"],
                "full_effect_rms": lolo["full_effect_rms"],
            }
        )
        observed_rms = vector_rms(observed)
        residual_rms = vector_rms(residual)
        pair_summary_rows.append(
            {
                "canonical_perturbation": key,
                "target_a": gene_a,
                "target_b": gene_b,
                "original_condition_labels": "|".join(originals[key]),
                "total_n": size_n[key],
                "minimum_lane_n": packed["minimum_lane_n"],
                "observed_rms": observed_rms,
                "additive_rmse": residual_rms,
                "delta_pearson": delta_pearson(observed, additive),
                "cosine": cosine_similarity(observed, additive),
                "residual_over_response_rms": zero_safe_ratio(residual_rms, observed_rms),
                "mean_same_sign_lane_fraction": float(np.mean(packed["same_sign_lane_fraction"])),
                "fraction_genes_majority_lane_sign": float(
                    np.mean(packed["same_sign_lane_fraction"] >= 0.5)
                ),
                "se_undefined_lane_n_lt_2": packed["se_undefined_lane_n_lt_2"],
                "lolo_max_rms_change": lolo["max_rms_change"],
            }
        )

    print("describe: writing long gene tables", flush=True)
    single_long = _long_frame(
        single_keys,
        genes,
        originals,
        single_mean,
        single_effect,
        single_se,
        single_lo,
        single_hi,
        single_pos,
        "positive_lane_fraction",
        single_min_n,
        [size_n[key] for key in single_keys],
        single_used,
        single_flag,
    )
    pair_long = _long_frame(
        pair_keys,
        genes,
        originals,
        pair_mean,
        pair_resid,
        pair_se,
        pair_lo,
        pair_hi,
        pair_sign,
        "same_sign_lane_fraction",
        pair_min_n,
        [size_n[key] for key in pair_keys],
        pair_used,
        pair_flag,
        extra_cols=(
            ("observed_delta", pair_obs),
            ("additive_delta", pair_add),
            ("residual", pair_resid),
        ),
    )
    single_long.to_csv(
        DESCRIBE_DIR / "single_effects_long.csv.gz",
        index=False,
        encoding="utf-8",
        compression="gzip",
    )
    pair_long.to_csv(
        DESCRIBE_DIR / "pair_interactions_long.csv.gz",
        index=False,
        encoding="utf-8",
        compression="gzip",
    )

    single_summary = pd.DataFrame(single_summary_rows).sort_values(
        "effect_rms",
        ascending=False,
    )
    pair_summary = pd.DataFrame(pair_summary_rows).sort_values(
        "additive_rmse",
        ascending=False,
    )
    single_summary.to_csv(DESCRIBE_DIR / "single_summaries.csv", index=False, encoding="utf-8")
    pair_summary.to_csv(DESCRIBE_DIR / "pair_summaries.csv", index=False, encoding="utf-8")
    pd.DataFrame(single_lolo + pair_lolo).to_csv(
        DESCRIBE_DIR / "leave_one_lane_out.csv",
        index=False,
        encoding="utf-8",
    )

    top_rows = []
    for i, key in enumerate(single_keys):
        order = np.argsort(-single_effect[i])
        for rank, idx in enumerate(order[:10], start=1):
            top_rows.append(
                {
                    "canonical_perturbation": key,
                    "direction": "up",
                    "rank": rank,
                    "gene_id": genes.loc[idx, "gene_id"],
                    "gene_name": genes.loc[idx, "gene_name"],
                    "effect": float(single_effect[i, idx]),
                }
            )
        for rank, idx in enumerate(order[-10:][::-1], start=1):
            top_rows.append(
                {
                    "canonical_perturbation": key,
                    "direction": "down",
                    "rank": rank,
                    "gene_id": genes.loc[idx, "gene_id"],
                    "gene_name": genes.loc[idx, "gene_name"],
                    "effect": float(single_effect[i, idx]),
                }
            )
    pd.DataFrame(top_rows).to_csv(
        DESCRIBE_DIR / "top_single_genes.csv",
        index=False,
        encoding="utf-8",
    )

    grouped = group_orientation_labels(
        [label for labels in originals.values() for label in labels]
    )
    dual = dual_orientation_single_keys(grouped)
    orientation_rows = []
    for key in dual:
        labels = grouped[key]
        pooled = effect_by_single[key]
        effects = []
        for label in labels:
            entry = original[label]
            contrasts = single_effect_lane(entry["mean"], ctrl["mean"])
            estimate = summarize_lane_contrasts(
                contrasts,
                [entry["n"], ctrl["n"]],
                [entry["sample_var"], ctrl["sample_var"]],
            )["estimate"]
            effects.append(estimate)
            orientation_rows.append(
                {
                    "canonical_perturbation": key,
                    "original_condition": label,
                    "rms_vs_pooled_single_effect": vector_rms(estimate - pooled),
                    "pooled_effect_rms": vector_rms(pooled),
                    "orientation_effect_rms": vector_rms(estimate),
                    "total_n_orientation": int(original[label]["n"].sum()),
                }
            )
        orientation_rows.append(
            {
                "canonical_perturbation": key,
                "original_condition": "BETWEEN_ORIENTATIONS",
                "rms_vs_pooled_single_effect": vector_rms(effects[0] - effects[1]),
                "pooled_effect_rms": vector_rms(pooled),
                "orientation_effect_rms": None,
                "total_n_orientation": None,
            }
        )
    orientation = pd.DataFrame(orientation_rows)
    orientation.to_csv(
        DESCRIBE_DIR / "orientation_sensitivity.csv",
        index=False,
        encoding="utf-8",
    )

    all_keys = ["ctrl"] + single_keys + pair_keys
    signature = np.stack(
        [np.zeros(n_genes, dtype=np.float64) if key == "ctrl" else condition_mean(canonical[key]) - ctrl_mean for key in all_keys]
    )
    svd = condition_signature_svd(signature, n_components=5)
    module_rows = []
    for pc in range(svd["n_components"]):
        loading = svd["loadings"][pc]
        order = np.argsort(-np.abs(loading))[:20]
        for rank, idx in enumerate(order, start=1):
            module_rows.append(
                {
                    "component": pc + 1,
                    "rank": rank,
                    "gene_id": genes.loc[idx, "gene_id"],
                    "gene_name": genes.loc[idx, "gene_name"],
                    "loading": float(loading[idx]),
                    "singular_value": float(svd["singular_values"][pc]),
                }
            )
    pd.DataFrame(module_rows).to_csv(
        DESCRIBE_DIR / "svd_modules.csv",
        index=False,
        encoding="utf-8",
    )
    pd.DataFrame(
        {
            "canonical_perturbation": all_keys,
            "n_targets": [0 if key == "ctrl" else 1 if "+" not in key else 2 for key in all_keys],
            "pc1": svd["scores"][:, 0],
            "pc2": svd["scores"][:, 1],
            "pc3": svd["scores"][:, 2],
        }
    ).to_csv(DESCRIBE_DIR / "condition_pca_scores.csv", index=False, encoding="utf-8")

    largest_pair = str(pair_summary.iloc[0]["canonical_perturbation"])
    largest_idx = pair_keys.index(largest_pair)
    plt = _try_pyplot()
    figure_status = {"backend": "matplotlib_agg", "generated": []}
    if plt is None:
        write_json(
            DESCRIBE_DIR / "figures_skipped.json",
            {"reason": "matplotlib_not_importable", "install_not_attempted": True},
        )
        figure_status["generated"] = []
        figure_status["skipped"] = True
    else:
        top_single = single_summary.head(20)["canonical_perturbation"].tolist()
        top_single_idx = [single_keys.index(key) for key in top_single]
        gene_score = np.mean(np.abs(single_effect[top_single_idx]), axis=0)
        top_genes = np.argsort(-gene_score)[:25]
        fig, ax = plt.subplots(figsize=(11, 7))
        panel = single_effect[np.ix_(top_single_idx, top_genes)]
        limit = float(np.max(np.abs(panel))) or 1.0
        im = ax.imshow(panel, aspect="auto", cmap="coolwarm", vmin=-limit, vmax=limit)
        ax.set_yticks(range(len(top_single)))
        ax.set_yticklabels(top_single, fontsize=7)
        ax.set_xticks(range(len(top_genes)))
        ax.set_xticklabels(genes.loc[top_genes, "gene_name"], rotation=90, fontsize=7)
        fig.colorbar(im, ax=ax, label="lane-balanced single effect (log1p X)")
        ax.set_title("Top single-effect signatures (descriptive, all cells)")
        _save_figure(
            fig,
            "single_effect_top_gene_heatmap",
            "Lane-balanced single-gene effects for the 20 singles with largest effect RMS, "
            "showing the 25 genes with largest mean absolute effect among those singles.",
            "Stored log1p expression difference versus all-data lane-balanced control.",
            "Full-data descriptive ranking, not a held-out test and not a causal claim. "
            "Eight lanes are technical GEM groups, not biological replicates.",
        )
        plt.close(fig)

        top_pair = pair_summary.head(15)["canonical_perturbation"].tolist()
        top_pair_idx = [pair_keys.index(key) for key in top_pair]
        inter_score = np.mean(np.abs(pair_resid[top_pair_idx]), axis=0)
        inter_genes = np.argsort(-inter_score)[:25]
        fig, ax = plt.subplots(figsize=(11, 6))
        panel = pair_resid[np.ix_(top_pair_idx, inter_genes)]
        limit = float(np.max(np.abs(panel))) or 1.0
        im = ax.imshow(panel, aspect="auto", cmap="coolwarm", vmin=-limit, vmax=limit)
        ax.set_yticks(range(len(top_pair)))
        ax.set_yticklabels(top_pair, fontsize=7)
        ax.set_xticks(range(len(inter_genes)))
        ax.set_xticklabels(genes.loc[inter_genes, "gene_name"], rotation=90, fontsize=7)
        fig.colorbar(im, ax=ax, label="lane-balanced interaction residual (log1p X)")
        ax.set_title("Strongest non-additive pair residuals (descriptive, all cells)")
        _save_figure(
            fig,
            "strongest_interaction_heatmap",
            "Lane-balanced pair residuals versus log-space additivity for the 15 pairs "
            "with largest residual RMS, showing 25 genes with largest mean |residual|.",
            "log1p residual = mean(pair)-mean(A)-mean(B)+mean(ctrl), then averaged across lanes.",
            "Additivity is on stored log1p space. This is not mechanistic synergy proof. "
            "Pair ranking used all cells and is exploratory.",
        )
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(6, 6))
        ax.scatter(pair_add[largest_idx], pair_obs[largest_idx], s=6, alpha=0.35, linewidths=0)
        limits = [
            min(pair_add[largest_idx].min(), pair_obs[largest_idx].min()),
            max(pair_add[largest_idx].max(), pair_obs[largest_idx].max()),
        ]
        ax.plot(limits, limits, color="black", linewidth=1)
        ax.set_xlabel("Additive delta (log1p X)")
        ax.set_ylabel("Observed pair delta (log1p X)")
        ax.set_title(f"Observed vs additive delta: {largest_pair}")
        _save_figure(
            fig,
            "observed_vs_additive_delta_scatter_largest_residual_pair",
            f"Per-gene observed pair delta versus additive prediction for {largest_pair}, "
            "the largest residual-RMS pair on the full-data descriptive ranking.",
            "Delta vs all-data lane-balanced control, stored log1p X.",
            "This pair was selected from the full-data exploratory ranking, not as a "
            "pre-registered held-out contrast. Do not read it as a paper-final story.",
        )
        plt.close(fig)

        fig = plot_pair_error_distribution_ranking(pair_summary, plt)
        _save_figure(
            fig,
            "pair_error_distribution_ranking",
            PAIR_RANKING_CAPTION,
            PAIR_RANKING_UNITS,
            PAIR_RANKING_CAVEAT,
        )
        plt.close(fig)

        test_metrics = metrics[metrics["split"] == "test"].copy()
        pivot = test_metrics.pivot(
            index="canonical_perturbation",
            columns="model",
            values="mse",
        )[list(BASELINE_ORDER)]
        fig, ax = plt.subplots(figsize=(12, 5))
        x = np.arange(len(pivot))
        width = 0.16
        for j, model in enumerate(BASELINE_ORDER):
            ax.bar(x + (j - 2) * width, pivot[model].to_numpy(), width=width, label=model)
        ax.set_xticks(x)
        ax.set_xticklabels(pivot.index, rotation=90, fontsize=7)
        ax.set_ylabel("Per-condition MSE")
        ax.set_title("Held-out test pair MSE by predeclared baseline")
        ax.legend(fontsize=6, loc="upper right")
        _save_figure(
            fig,
            "heldout_baseline_paired_condition_error",
            "Paired per-condition MSE on the 26 B-protocol test pairs for the five "
            "predeclared baselines. Ridge alpha was frozen from validation only.",
            "Macro-equal-condition MSE on 5,045-gene deltas versus train-control mu0.",
            "Test was evaluated once. This is predictive combination generalization, "
            "not a causal interaction identified by the split.",
        )
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(6, 5))
        n_targets_color = [0 if key == "ctrl" else 1 if "+" not in key else 2 for key in all_keys]
        scatter = ax.scatter(
            svd["scores"][:, 0],
            svd["scores"][:, 1],
            c=n_targets_color,
            cmap="viridis",
            s=18,
        )
        fig.colorbar(scatter, ax=ax, label="n targets (0=ctrl, 1=single, 2=pair)")
        ax.set_xlabel("PC1 (all-data SVD)")
        ax.set_ylabel("PC2 (all-data SVD)")
        ax.set_title("237 condition signatures (exploratory SVD)")
        _save_figure(
            fig,
            "pca_237_condition_signatures_svd",
            "First two principal coordinates from a column-centered SVD of 237 "
            "lane-balanced condition deltas, including control as the origin.",
            "Scores are in log1p-delta space. Not used as training features.",
            "Exploratory all-data geometry. Not a UMAP. Full-cohort HVG caveat applies.",
        )
        plt.close(fig)
        figure_status["generated"] = [
            "single_effect_top_gene_heatmap",
            "strongest_interaction_heatmap",
            "observed_vs_additive_delta_scatter_largest_residual_pair",
            "pair_error_distribution_ranking",
            "heldout_baseline_paired_condition_error",
            "pca_237_condition_signatures_svd",
        ]
        figure_status["skipped"] = False

    write_json(DESCRIBE_DIR / "figure_status.json", figure_status)
    write_json(
        DESCRIBE_DIR / "describe_meta.json",
        {
            "used_in_fit": False,
            "chosen_alpha_frozen": chosen["alpha"],
            "n_single_tables": 105,
            "n_pair_tables": 131,
            "n_genes": n_genes,
            "largest_residual_pair_full_data": largest_pair,
            "n_dual_orientation_singles": len(dual),
            "go_enrichment": "skipped",
            "significance": SIGNIFICANCE_NOT_ESTIMABLE,
        },
    )
    write_modeling_report(
        protocol,
        prepare_meta,
        chosen,
        fit_summaries,
        bootstrap,
        examples,
        sizes,
        single_summary,
        pair_summary,
        orientation,
        curve,
        metrics,
        figure_status,
    )
    print("describe: wrote report " + str(REPORT_PATH), flush=True)
    return {"largest_residual_pair_full_data": largest_pair, "chosen_alpha_frozen": chosen["alpha"]}


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
    body = []
    for row in rows:
        body.append("| " + " | ".join(_fmt(row.get(col, "")).replace("|", "; ").replace("\n", " ") for col in columns) + " |")
    return "\n".join([header, sep, *body])


def write_modeling_report(
    protocol: dict,
    prepare_meta: dict,
    chosen: dict,
    summaries: dict,
    bootstrap: dict,
    examples: pd.DataFrame,
    sizes: pd.DataFrame,
    single_summary: pd.DataFrame,
    pair_summary: pd.DataFrame,
    orientation: pd.DataFrame,
    curve: pd.DataFrame,
    metrics: pd.DataFrame,
    figure_status: dict,
) -> None:
    test = summaries["test"]
    val = summaries["validation"]
    top_singles = single_summary.head(5)
    top_pairs = pair_summary.head(5)
    example_single = single_summary.iloc[0]
    example_pair = pair_summary.iloc[0]
    ori_between = orientation.loc[orientation["original_condition"] == "BETWEEN_ORIENTATIONS"]
    max_ori = float(ori_between["rms_vs_pooled_single_effect"].max()) if len(ori_between) else float("nan")
    test_ridge = test["ridge_multihot_main_effects"]
    test_add = test["additive"]
    recommend = _recommend_next(test)

    example_rows = []
    for row in examples.itertuples(index=False):
        example_rows.append(
            {
                "cell_barcode": row.cell_barcode,
                "condition": row.condition,
                "canonical": row.canonical_perturbation,
                "lane": row.barcode_group,
                "B_split": getattr(row, protocol["split"]["column"]),
            }
        )

    largest_sizes = sizes.sort_values("total_n", ascending=False).head(8)
    size_rows = [
        {
            "canonical_perturbation": row.canonical_perturbation,
            "total_n": int(row.total_n),
            "n_targets": int(row.n_targets),
            "labels": row.original_condition_labels,
        }
        for row in largest_sizes.itertuples(index=False)
    ]
    single_rows = [
        {
            "perturbation": row.canonical_perturbation,
            "effect_rms": row.effect_rms,
            "total_n": int(row.total_n),
            "lolo_max_rms_change": row.lolo_max_rms_change,
            "top10_up": row.top10_up_genes,
            "top10_down": row.top10_down_genes,
        }
        for row in top_singles.itertuples(index=False)
    ]
    pair_rows = [
        {
            "pair": row.canonical_perturbation,
            "observed_rms": row.observed_rms,
            "additive_rmse": row.additive_rmse,
            "cosine": row.cosine,
            "delta_pearson": row.delta_pearson,
            "resid/response": row.residual_over_response_rms,
            "sign_support": row.mean_same_sign_lane_fraction,
        }
        for row in top_pairs.itertuples(index=False)
    ]
    curve_rows = [
        {"alpha": row.alpha, "validation_macro_mse": row.validation_macro_mse}
        for row in curve.itertuples(index=False)
    ]
    model_rows = []
    for model in BASELINE_ORDER:
        model_rows.append(
            {
                "model": model,
                "val_macro_mse": val[model]["macro_mse"],
                "test_macro_mse": test[model]["macro_mse"],
                "test_macro_mae": test[model]["macro_mae"],
                "test_delta_pearson": test[model]["macro_delta_pearson"],
                "test_cosine": test[model]["macro_cosine"],
            }
        )
    example_table = _md_table(
        example_rows,
        ["cell_barcode", "condition", "canonical", "lane", "B_split"],
    )
    size_table = _md_table(
        size_rows,
        ["canonical_perturbation", "total_n", "n_targets", "labels"],
    )
    single_table = _md_table(
        single_rows,
        ["perturbation", "effect_rms", "total_n", "lolo_max_rms_change", "top10_up", "top10_down"],
    )
    pair_table = _md_table(
        pair_rows,
        ["pair", "observed_rms", "additive_rmse", "cosine", "delta_pearson", "resid/response", "sign_support"],
    )
    curve_table = _md_table(curve_rows, ["alpha", "validation_macro_mse"])
    score_table = _md_table(
        model_rows,
        ["model", "val_macro_mse", "test_macro_mse", "test_macro_mae", "test_delta_pearson", "test_cosine"],
    )

    text = f"""# Modeling report — Norman / GEARS-processed Perturb-seq

Generated by `scripts/modeling_analysis.py` from files under `reports/modeling/`.
Numbers below are copied from those artifacts. This report does not invent biological examples.

**Citations.** Norman et al., *Science* 2019, doi:[10.1126/science.aax4438](https://doi.org/10.1126/science.aax4438). Protocol decisions and baseline roles: [METHODS_DECISION_BRIEF.md](METHODS_DECISION_BRIEF.md). Dataset identity and split provenance: [DATA_PROVENANCE.md](DATA_PROVENANCE.md). Predeclared protocol: [../config/modeling_protocol.yaml](../config/modeling_protocol.yaml).

## Distinctions that apply to every section

| Class | What it is | What it is not |
| --- | --- | --- |
| Observed / descriptive | All-cell lane-balanced tables written after fit | A held-out score; a reason to change alpha |
| Held-out test | 26 B-protocol pairs, evaluated once after validation freeze | A biological confidence interval; a causal interaction test |
| Unknown | Orientation interchangeability, mechanistic synergy, donor-level uncertainty | Resolved by eight technical lanes |

Uncertainty is a **conditional cell-sampling** SE under disjoint lanes and conditionally independent captured cells. The significance field is `{SIGNIFICANCE_NOT_ESTIMABLE}`. Eight GEM lanes are not eight biological replicates. Rankings of observed effects are descriptive. Uncertainty is not proof of causality. Additivity is on **stored log1p X**, not a mechanistic synergy proof. Construct orientations that share a target set were pooled; that is a caveat, not an experimental equivalence proof.

## 1. Dataset structure

What can be modeled is the stored 91,205 × 5,045 log1p matrix `X` (CSR; counts payload was not loaded). Feature order is immutable `gene_id`. Stored `cell_type=A549` remains a dummy; the paper/GEO identity is K562 SunTag CRISPRa.

| Quantity | Value |
| --- | ---: |
| Cells | {prepare_meta["n_obs"]} |
| Measured genes | {prepare_meta["n_vars"]} |
| Original condition labels | {prepare_meta["n_original_condition_labels"]} |
| Canonical target sets including control | {prepare_meta["n_canonical_target_sets_including_control"]} |
| Control / single / pair cells | {prepare_meta["n_control_cells"]} / {prepare_meta["n_single_cells"]} / {prepare_meta["n_double_cells"]} |
| Dual-orientation single target sets | {len(prepare_meta["dual_orientation_singles"])} |
| `X` nonzeros | {prepare_meta["x_nnz"]} |
| Counts payload loaded | {prepare_meta["counts_payload_loaded"]} |
| `uns` read | {prepare_meta["uns_read"]} |

Available representations: `layers[counts]` contains integer-valued raw UMI counts for the retained features; `X` contains log1p-transformed per-cell normalized values. There is no separately stored unlogged normalized layer and no original PCA/latent embedding. The PCA below is newly computed on condition signatures, not individual cells. Original obs fields are condition, cell_type, dose_val, control, and condition_name. Barcode suffixes encode technical lanes. Guide identity and guide-assignment QC are recovered in [verified_cell_metadata.csv](verified_cell_metadata.csv), rather than stored in obs. No donor or independent biological-replicate ID is available. Controls carry non-targeting constructs, rather than representing untreated tracked twins.

The 5,045 genes are a **full-cohort prepared HVG-plus-targets set**. They are usable retrospectively with that caveat and were not re-selected here. Stored `uns` DE lists and external `bulk_fitness` were not used.

First stored cells (not a matrix dump):

{example_table}

Largest canonical conditions by cell count:

{size_table}

Full sizes: [modeling/prepare/condition_sizes.csv](modeling/prepare/condition_sizes.csv).

## 2. Target definition

Predictive target: the 5,045-gene **lane-balanced condition mean** of stored log1p `X`, as a delta versus **training-split controls only** (`mu0`). Validation/test controls are diagnostic and were not model inputs or the reference.

Descriptive target (computed after fit): the same lane-balanced mean, but using **all** cells of a condition and all-data controls. That analysis cannot change the frozen ridge alpha `{_fmt(chosen["alpha"])}`.

Single-gene effect in lane ℓ: `mean(single,ℓ) − mean(ctrl,ℓ)`, then the eight-lane average. Pair interaction in lane ℓ: `mean(pair,ℓ) − mean(A,ℓ) − mean(B,ℓ) + mean(ctrl,ℓ)`, then the eight-lane average. SE uses sample variance / n per lane mean, sums contrast variances, and divides by 64. If any component lane has n<2, SE is NA.

The response vector retains one effect per measured gene and makes single/pair contrasts interpretable. It averages cell-wise logged expression; it is **not** log of summed UMI counts or count-pseudobulk differential expression. Count-pseudobulk can support a different estimand, but this experiment lacks independent biological replicates for replicate-based DE. Cell-level distributions can reveal heterogeneity hidden by these means, at greater estimation and evaluation cost; unpaired cells still cannot identify individual counterfactual effects. No distribution model was fit.

## 3. Individual perturbation effects

All-data descriptive. Full long table (105 × 5,045): [modeling/describe/single_effects_long.csv.gz](modeling/describe/single_effects_long.csv.gz). Summaries: [modeling/describe/single_summaries.csv](modeling/describe/single_summaries.csv). Top 10 up/down genes per single: [modeling/describe/top_single_genes.csv](modeling/describe/top_single_genes.csv).

Largest effect-RMS singles:

{single_table}

Example from the tables: perturbing **{example_single.canonical_perturbation}** has effect RMS {_fmt(example_single.effect_rms)} on stored log1p X (n={int(example_single.total_n)} cells). Strongest up genes in that row: {example_single.top10_up_genes}. Strongest down genes: {example_single.top10_down_genes}. Leave-one-lane-out max RMS change is {_fmt(example_single.lolo_max_rms_change)}. This is an observed signature, not a causal proof.

## 4. Combination analysis

All-data descriptive. Long table (131 × 5,045): [modeling/describe/pair_interactions_long.csv.gz](modeling/describe/pair_interactions_long.csv.gz). Pair summaries: [modeling/describe/pair_summaries.csv](modeling/describe/pair_summaries.csv).

Largest additive residual RMS pairs:

{pair_table}

Example from the tables: **{example_pair.canonical_perturbation}** has observed RMS {_fmt(example_pair.observed_rms)}, additive RMSE {_fmt(example_pair.additive_rmse)}, cosine {_fmt(example_pair.cosine)}, delta Pearson {_fmt(example_pair.delta_pearson)}, and residual/response RMS {_fmt(example_pair.residual_over_response_rms)}. These ranks used all cells.

## 5. Interaction effects

The interaction residual is the log-space quantity above. Strongest non-additive pairs are the residual-RMS ranking in the pair table. Formal GO enrichment was **skipped** (no validated broad measured-gene mapping was used here). Descriptive expression modules from SVD loadings are in [modeling/describe/svd_modules.csv](modeling/describe/svd_modules.csv). They are not training features.

## 6. Baseline modeling results

Split: existing `B_heldout_combination` seed {protocol["split"]["seed"]} exactly (105 singles + 79 pairs train; 26 / 26 pairs val/test). All centroids and `mu0` used training cells only. Ridge features are 105-gene multi-hot main effects, no pair one-hots. Alpha grid {protocol["ridge"]["alphas"]}; one alpha selected by validation macro MSE, then frozen. Train+val was not refit. Test was evaluated once.

Chosen alpha: **{_fmt(chosen["alpha"])}**.

Validation curve:

{curve_table}

Scorecard (primary = condition-macro MSE; MAE secondary; Pearson NA for constant vectors, cosine NA for zero-norm vectors):

{score_table}

Matching-single and additive predictions differ by a positive factor of two in delta space; their Pearson/cosine scores are therefore identical by construction. MSE distinguishes their magnitude errors.

Test skills for ridge: vs perturbed mean {_fmt(test_ridge.get("skill_vs_train_condition_balanced_perturbed_mean"))}; vs additive {_fmt(test_ridge.get("skill_vs_additive"))}. Additive vs perturbed mean {_fmt(test_add.get("skill_vs_train_condition_balanced_perturbed_mean"))}.

Paired bootstrap of test-condition MSE(ridge) − MSE(additive), seed {bootstrap["seed"]}, {bootstrap["n_samples"]} samples, **conditional on the frozen model and train reference**: mean difference {_fmt(bootstrap["mean_difference"])}; 2.5/97.5 percentiles {_fmt(bootstrap["percentile_2_5"])} / {_fmt(bootstrap["percentile_97_5"])}. This is not a biological CI. Pairs share target genes and single-effect reference estimates, so resampling conditions omits training/reference uncertainty and must not be read as population-level significance.

Val/test control profiles vs train `mu0` were recorded only as diagnostics in [modeling/fit/chosen_alpha.json](modeling/fit/chosen_alpha.json).

Per-condition errors: [modeling/fit/per_condition_metrics.csv](modeling/fit/per_condition_metrics.csv). Predictions: [modeling/fit/predictions.npz](modeling/fit/predictions.npz).

## 7. Evaluation design

Protocol B asks whether a map from seen singles and some seen pairs can forecast an unseen pair's mean log1p profile. The split does **not** identify a paired treatment effect. Random cell split A was not used. No test-expression graphs, no `uns` DE, no `bulk_fitness`. Ridge cannot represent pair-specific interactions that do not project onto main-effect multi-hots.

## 8. Recommended next model

No nonlinear model was trained (only 79 training pair conditions; linear/additive were required first).

{recommend}

## Orientation sensitivity and technical stability

{len(prepare_meta["dual_orientation_singles"])} single target sets have both stored orientations. Maximum RMS difference between orientations (descriptive, not a new biology claim): {_fmt(max_ori)}. Table: [modeling/describe/orientation_sensitivity.csv](modeling/describe/orientation_sensitivity.csv).

Leave-one-lane-out RMS ranges: [modeling/describe/leave_one_lane_out.csv](modeling/describe/leave_one_lane_out.csv). These assess technical-lane stability, not biological replication.

## Figures

Diagnostic plots used matplotlib Agg. They are not a frozen paper story. Each file has a sidecar caption with units and caveats.

Generated (click through; captions sit beside each PNG):

- [single-effect heatmap](modeling/describe/figures/single_effect_top_gene_heatmap.png)
- [strongest-interaction heatmap](modeling/describe/figures/strongest_interaction_heatmap.png) (zero-centered)
- [observed vs additive scatter](modeling/describe/figures/observed_vs_additive_delta_scatter_largest_residual_pair.png)
- [pair residual ranking](modeling/describe/figures/pair_error_distribution_ranking.png)
- [held-out baseline errors](modeling/describe/figures/heldout_baseline_paired_condition_error.png)
- [condition-signature PCA](modeling/describe/figures/pca_237_condition_signatures_svd.png)

Audited biological summary from diagnostic tables: [modeling/AUDITED_SUMMARY.md](modeling/AUDITED_SUMMARY.md). That file is written by `scripts/finalize_modeling_report.py` and is not used in fit.

## Repeatability

Deterministic given the local H5AD, [proposed_split_assignments.csv](proposed_split_assignments.csv), this protocol file, and the stated seeds (split 20260911 already in the manifest; bootstrap {BOOTSTRAP_SEED}). Re-running `prepare`/`fit`/`describe` in that order should regenerate the same chosen alpha and the same train predictions; test outcomes cannot change them. `describe` reads the frozen alpha and does not refit.

## Unknown

Whether `GENE+ctrl` and `ctrl+GENE` are experimentally interchangeable; whether non-additivity is mechanistic synergy; donor-level uncertainty; GEO raw UMI reconstruction; `SLC38A2` guide sequence; RHOXF2 / RHOXF2B / RHOXF2BB molecular identity.
"""
    REPORT_PATH.write_text(text, encoding="utf-8")


def _recommend_next(test: dict) -> str:
    ridge = test["ridge_multihot_main_effects"]["macro_mse"]
    additive = test["additive"]["macro_mse"]
    return (
        f"Keep ridge as the current predictive baseline (MSE {_fmt(ridge)} versus "
        f"additive {_fmt(additive)}). Its main-effect encoding cannot learn arbitrary "
        "pair-specific interactions: improvement can reflect shrinkage, recalibration, "
        "and shared response structure. These results alone do not justify a nonlinear "
        "or deep model. First examine construct sensitivity and shared-shift diagnostics "
        "in [MODELING_DIAGNOSTICS.md](MODELING_DIAGNOSTICS.md). A future interaction "
        "model needs a prespecified low-dimensional representation and an honest new "
        "evaluation. The 26 test pairs have now been inspected; do not reuse them for "
        "tuning. Use nested condition-level resampling as an explicitly retrospective "
        "benchmark, or obtain new independent held-out experiments for stronger claims."
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage",
        required=False,
        choices=("prepare", "fit", "describe", "all"),
        help="prepare, fit, describe, or all (fit before describe)",
    )
    parser.add_argument(
        "--plots-only",
        nargs="?",
        const="pair_error_distribution_ranking",
        default=None,
        help="Refresh one saved figure from pair_summaries; default stem is pair_error_distribution_ranking",
    )
    parser.add_argument(
        "--protocol",
        type=Path,
        default=DEFAULT_PROTOCOL,
        help="Path to modeling_protocol.yaml",
    )
    args = parser.parse_args()
    if args.plots_only:
        if args.plots_only != "pair_error_distribution_ranking":
            raise SystemExit("only pair_error_distribution_ranking can be refreshed without describe")
        path = refresh_pair_error_distribution_ranking()
        print("plots-only: wrote " + str(path), flush=True)
        return 0
    if not args.stage:
        raise SystemExit("specify --stage or --plots-only")
    protocol = load_protocol(args.protocol)
    stages = ("prepare", "fit", "describe") if args.stage == "all" else (args.stage,)
    if args.stage == "describe" and not (FIT_DIR / "chosen_alpha.json").is_file():
        raise SystemExit("describe requires a completed fit; it cannot select alpha")
    if "describe" in stages and "fit" not in stages and args.stage == "all":
        raise SystemExit("internal stage order error")
    for stage in stages:
        print(f"=== {stage} ===", flush=True)
        if stage == "prepare":
            stage_prepare(protocol)
        elif stage == "fit":
            stage_fit(protocol)
        else:
            stage_describe(protocol)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
