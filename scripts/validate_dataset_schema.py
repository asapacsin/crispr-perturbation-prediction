"""Read-only validation of config/dataset_schema.yaml against local artifacts.

Does not load the dense counts layer. Does not train models or mutate the H5AD.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import h5py
import pandas as pd

try:
    import yaml
except ImportError as exc:  # pragma: no cover - environment check
    raise SystemExit(
        "ERROR: PyYAML is required to read config/dataset_schema.yaml. "
        "The project brief said it is already installed; this script will not pip-install it. "
        f"Original import error: {exc}"
    ) from exc

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

from dataset_schema_lib import (
    BRIDGE_SPLIT,
    barcode_split_errors,
    canonicalize_condition,
    combination_split_errors,
    gene_partition_errors,
    n_targets,
    target_gene_split_errors,
)

DEFAULT_SCHEMA = ROOT / "config" / "dataset_schema.yaml"


class SchemaValidationError(Exception):
    pass


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _as_strings(values) -> list[str]:
    out = []
    for value in values:
        if isinstance(value, bytes):
            out.append(value.decode("utf-8"))
        else:
            out.append(str(value))
    return out


def _map_categorical_codes(categories: list[str], codes, *, name: str) -> list[str]:
    n_categories = len(categories)
    mapped: list[str] = []
    for raw in codes:
        code = int(raw)
        if code == -1:
            raise SchemaValidationError(
                f"null categorical code (-1) in {name}; required annotations cannot be missing"
            )
        if code < 0 or code >= n_categories:
            raise SchemaValidationError(
                f"categorical code {code} in {name} is out of bounds for {n_categories} categories"
            )
        mapped.append(categories[code])
    return mapped


def _categories_from_legacy_or_group(node, *, name: str) -> list[str]:
    if isinstance(node, h5py.Dataset):
        return _as_strings(node[()])
    if isinstance(node, h5py.Group) and "categories" in node:
        return _as_strings(node["categories"][()])
    raise SchemaValidationError(
        f"cannot read categories for {name}; "
        f"type={type(node).__name__}"
    )


def _read_anndata_vector(container: h5py.Group, name: str) -> list[str]:
    if name not in container:
        raise SchemaValidationError(f"{container.name} is missing {name}")
    node = container[name]
    path = f"{container.name}/{name}"
    if isinstance(node, h5py.Dataset):
        values = node[()]
        categories_root = container.get("__categories")
        if categories_root is not None and name in categories_root:
            categories = _categories_from_legacy_or_group(
                categories_root[name],
                name=f"{categories_root.name}/{name}",
            )
            return _map_categorical_codes(categories, values, name=path)
        return _as_strings(values)
    if isinstance(node, h5py.Group):
        if "categories" in node and "codes" in node:
            categories = _as_strings(node["categories"][()])
            return _map_categorical_codes(categories, node["codes"][()], name=path)
    raise SchemaValidationError(
        f"cannot decode AnnData vector {path}; "
        f"type={type(node).__name__} keys={list(node.keys()) if isinstance(node, h5py.Group) else None}"
    )


def _read_numeric(container: h5py.Group, name: str):
    if name not in container:
        raise SchemaValidationError(f"{container.name} is missing {name}")
    node = container[name]
    if isinstance(node, h5py.Dataset):
        return node[()]
    if isinstance(node, h5py.Group) and "codes" in node:
        return node["codes"][()]
    raise SchemaValidationError(f"cannot decode numeric vector at {node.name}")


def matrix_shape(node) -> tuple[int, int]:
    if isinstance(node, h5py.Dataset):
        return tuple(int(x) for x in node.shape)
    attrs = dict(node.attrs)
    if "shape" in attrs:
        return tuple(int(x) for x in attrs["shape"])
    if "h5sparse_shape" in attrs:
        return tuple(int(x) for x in attrs["h5sparse_shape"])
    if "shape" in node:
        return tuple(int(x) for x in node["shape"][()])
    raise SchemaValidationError(f"no shape recorded for {node.name}")


def load_obs_frame(h5: h5py.File) -> pd.DataFrame:
    obs = h5["obs"]
    index_name = "_index" if "_index" in obs else None
    if index_name is None and "cell_barcode" in obs:
        index_name = "cell_barcode"
    if index_name is None:
        raise SchemaValidationError(
            f"HDF5 obs is missing _index; present={list(obs.keys())}"
        )
    barcodes = _read_anndata_vector(obs, index_name)
    frame = pd.DataFrame({"cell_barcode": barcodes})
    required = ("condition", "cell_type", "dose_val", "control", "condition_name")
    missing = [name for name in required if name not in obs]
    if missing:
        raise SchemaValidationError(
            f"HDF5 obs missing columns {missing}; present={list(obs.keys())}"
        )
    frame["condition"] = _read_anndata_vector(obs, "condition")
    frame["cell_type"] = _read_anndata_vector(obs, "cell_type")
    frame["dose_val"] = _read_anndata_vector(obs, "dose_val")
    frame["condition_name"] = _read_anndata_vector(obs, "condition_name")
    frame["control"] = [int(v) for v in _read_numeric(obs, "control")]
    return frame


def load_var_frame(h5: h5py.File) -> pd.DataFrame:
    var = h5["var"]
    index_name = "_index" if "_index" in var else None
    if index_name is None and "gene_id" in var:
        index_name = "gene_id"
    if index_name is None:
        raise SchemaValidationError(
            f"HDF5 var is missing _index; present={list(var.keys())}"
        )
    if "gene_name" not in var:
        raise SchemaValidationError("HDF5 var is missing gene_name")
    return pd.DataFrame(
        {
            "gene_id": _read_anndata_vector(var, index_name),
            "gene_name": _read_anndata_vector(var, "gene_name"),
        }
    )


def add_error(errors: list[str], message: str) -> None:
    errors.append(message)


def validate(schema_path: Path = DEFAULT_SCHEMA) -> list[str]:
    errors: list[str] = []
    schema = yaml.safe_load(schema_path.read_text(encoding="utf-8"))
    if schema.get("final_protocol_approved") is not False:
        add_error(errors, "schema.final_protocol_approved must be false until the user approves a protocol")
    if schema.get("training_allowed") is not False:
        add_error(errors, "schema.training_allowed must be false; this schema does not authorize training")

    identity = schema["identity"]
    h5ad_path = ROOT / identity["local_h5ad"]
    if not h5ad_path.is_file():
        add_error(errors, f"H5AD not found: {h5ad_path}")
        return errors

    actual_sha = file_sha256(h5ad_path)
    expected_sha = identity["sha256"]
    if actual_sha != expected_sha:
        add_error(
            errors,
            f"SHA256 mismatch for {h5ad_path}: expected {expected_sha}, got {actual_sha}",
        )

    with h5py.File(h5ad_path, "r") as h5:
        expected_keys = set(identity["hdf5_root_keys"])
        actual_keys = set(h5.keys())
        if actual_keys != expected_keys:
            add_error(
                errors,
                f"HDF5 root keys {sorted(actual_keys)} != schema {sorted(expected_keys)}",
            )
        if "layers" not in h5:
            add_error(errors, "HDF5 has no layers group")
        else:
            layer_keys = list(h5["layers"].keys())
            if layer_keys != ["counts"]:
                add_error(
                    errors,
                    f"stored layers are {layer_keys}; schema allows only ['counts']. "
                    "Do not load the counts payload during this check.",
                )
            else:
                counts_shape = matrix_shape(h5["layers"]["counts"])
                expected_shape = tuple(identity["shape"])
                if counts_shape != expected_shape:
                    add_error(
                        errors,
                        f"counts shape {counts_shape} != schema {expected_shape}",
                    )
        if "X" not in h5:
            add_error(errors, "HDF5 is missing X")
        else:
            x_shape = matrix_shape(h5["X"])
            expected_shape = tuple(identity["shape"])
            if x_shape != expected_shape:
                add_error(errors, f"X shape {x_shape} != schema {expected_shape}")
        if "raw" in h5:
            add_error(errors, "HDF5 unexpectedly contains raw; schema says none is stored")
        obs = load_obs_frame(h5)
        var = load_var_frame(h5)

    if len(obs) != identity["n_obs"]:
        add_error(errors, f"n_obs {len(obs)} != schema {identity['n_obs']}")
    if len(var) != identity["n_vars"]:
        add_error(errors, f"n_vars {len(var)} != schema {identity['n_vars']}")
    if not obs["cell_barcode"].is_unique:
        add_error(errors, "cell_barcode values are not unique")
    if not var["gene_id"].is_unique or not var["gene_name"].is_unique:
        add_error(errors, "gene_id or gene_name is not unique")

    try:
        obs["canonical_perturbation"] = obs["condition"].map(canonicalize_condition)
        obs["n_targets"] = obs["condition"].map(n_targets)
    except Exception as exc:
        add_error(errors, f"condition parsing failed: {exc}")
        return errors

    expected_targets = list(schema["targets"]["names"])
    if len(expected_targets) != schema["targets"]["n"]:
        add_error(
            errors,
            f"schema targets.n is {schema['targets']['n']} but names has {len(expected_targets)}",
        )
    if len(set(expected_targets)) != len(expected_targets):
        add_error(errors, "schema target names are not unique")
    observed_targets = sorted(
        {gene for label in obs.loc[obs["n_targets"] > 0, "canonical_perturbation"] for gene in label.split("+")}
    )
    if observed_targets != expected_targets:
        missing = sorted(set(expected_targets) - set(observed_targets))
        extra = sorted(set(observed_targets) - set(expected_targets))
        add_error(
            errors,
            f"parsed target names != schema 105-name list; missing={missing[:12]} extra={extra[:12]}",
        )

    var_names = set(var["gene_name"])
    unmatched = [name for name in expected_targets if name not in var_names]
    if unmatched:
        add_error(errors, f"schema targets missing from var.gene_name: {unmatched}")

    id_map = {row["name"]: row["gene_id"] for row in schema["targets"]["gene_ids"]}
    for name, gene_id in id_map.items():
        hit = var.loc[var["gene_name"] == name, "gene_id"]
        if hit.empty or hit.iloc[0] != gene_id:
            add_error(errors, f"target {name} gene_id {gene_id} does not match HDF5 var")

    control = schema["control"]
    control_mask = obs["condition"] == control["condition_value"]
    if int(control_mask.sum()) != control["n_cells"]:
        add_error(
            errors,
            f"control condition count {int(control_mask.sum())} != schema {control['n_cells']}",
        )
    if not (control_mask == (obs["control"] == control["control_value"])).all():
        add_error(errors, "obs.control does not agree with condition == ctrl for every row")

    counts = schema["perturbation_counts"]
    by_n = obs["n_targets"].value_counts().to_dict()
    if by_n.get(0, 0) != counts["control_cells"]:
        add_error(errors, f"control cells {by_n.get(0, 0)} != schema {counts['control_cells']}")
    if by_n.get(1, 0) != counts["single_cells"]:
        add_error(errors, f"single cells {by_n.get(1, 0)} != schema {counts['single_cells']}")
    if by_n.get(2, 0) != counts["double_cells"]:
        add_error(errors, f"double cells {by_n.get(2, 0)} != schema {counts['double_cells']}")
    if obs["condition"].nunique() != counts["original_condition_labels"]:
        add_error(
            errors,
            f"original condition labels {obs['condition'].nunique()} != "
            f"{counts['original_condition_labels']}",
        )
    if obs["canonical_perturbation"].nunique() != counts["target_sets_including_control"]:
        add_error(
            errors,
            f"target sets {obs['canonical_perturbation'].nunique()} != "
            f"{counts['target_sets_including_control']}",
        )

    obs["barcode_group"] = obs["cell_barcode"].str.rsplit("-", n=1).str[-1].astype(int)
    expected_groups = {int(row["barcode_group"]): row for row in schema["groups"]["lanes"]}
    actual_groups = obs["barcode_group"].value_counts().to_dict()
    if set(actual_groups) != set(expected_groups):
        add_error(
            errors,
            f"barcode groups {sorted(actual_groups)} != schema {sorted(expected_groups)}",
        )
    for group, expected in expected_groups.items():
        actual = int(actual_groups.get(group, 0))
        if actual != expected["cells"]:
            add_error(errors, f"group {group} has {actual} cells, schema says {expected['cells']}")

    metadata_path = ROOT / schema["source_provenance"]["verified_cell_metadata"]
    if metadata_path.is_file():
        meta = pd.read_csv(metadata_path)
        if len(meta) != len(obs):
            add_error(
                errors,
                f"verified_cell_metadata has {len(meta)} rows, HDF5 obs has {len(obs)}",
            )
        merged = obs.merge(meta, on="cell_barcode", how="left", suffixes=("", "_src"))
        if merged["gemgroup"].isna().any():
            add_error(errors, "one or more HDF5 barcodes are missing from verified_cell_metadata")
        if not (merged["barcode_group"] == merged["gemgroup"]).all():
            add_error(errors, "barcode suffix does not match source gemgroup for every joined cell")
        coverage = merged["good_coverage"]
        if coverage.dtype == object:
            coverage = coverage.astype(str).str.lower().isin(["true", "1", "yes"])
        if not bool(coverage.all()) or not (pd.to_numeric(merged["number_of_cells"]) == 1).all():
            add_error(errors, "joined source rows are not all good_coverage=true and number_of_cells=1")
    else:
        add_error(errors, f"missing verified cell metadata: {metadata_path}")

    proposal = schema["proposal_only"]
    split_path = ROOT / proposal["split_manifest"]
    evidence_path = ROOT / proposal["split_evidence"]
    if not split_path.is_file():
        add_error(errors, f"missing proposed split manifest: {split_path}")
        return errors
    splits = pd.read_csv(split_path)
    if "cell_barcode" not in splits.columns:
        add_error(errors, "proposed_split_assignments.csv must contain cell_barcode")
        return errors
    if set(splits["cell_barcode"]) != set(obs["cell_barcode"]):
        add_error(errors, "proposed split barcodes do not exactly match HDF5 barcodes")
    split_obs = obs.merge(splits, on="cell_barcode", how="left", suffixes=("", "_split"))

    protocol_columns = {
        "A_random_cell": proposal["protocols"]["A_random_cell"]["column"],
        "B_heldout_combination": proposal["protocols"]["B_heldout_combination"]["column"],
        "C_heldout_target_gene": proposal["protocols"]["C_heldout_target_gene"]["column"],
        "D_heldout_barcode_group": proposal["protocols"]["D_heldout_barcode_group"]["column"],
    }
    for protocol, column in protocol_columns.items():
        if column not in split_obs.columns:
            add_error(errors, f"{protocol} column {column} missing from split manifest")
            continue
        assigned = list(zip(split_obs["cell_barcode"], split_obs[column].astype(str)))
        allowed = ("train", "validation", "test")
        if protocol == "C_heldout_target_gene":
            allowed = allowed + (BRIDGE_SPLIT,)
        for message in barcode_split_errors(assigned, allowed):
            add_error(errors, f"{protocol}: {message}")
        sizes = (
            split_obs.groupby(column)
            .size()
            .to_dict()
        )
        expected_sizes = proposal["protocols"][protocol]["sizes"]
        for split_name, expected in expected_sizes.items():
            actual = int(sizes.get(split_name, 0))
            if actual != expected["cells"]:
                add_error(
                    errors,
                    f"{protocol} {split_name} has {actual} cells, schema says {expected['cells']}",
                )

    b_rows = [
        (row.condition, row.B_heldout_combination)
        for row in split_obs.itertuples(index=False)
        if row.B_heldout_combination in {"train", "validation", "test"}
    ]
    for message in combination_split_errors(b_rows):
        add_error(errors, f"B_heldout_combination: {message}")

    c_schema = proposal["protocols"]["C_heldout_target_gene"]
    train_genes = set(c_schema["train_genes"])
    val_genes = set(c_schema["validation_genes"])
    test_genes = set(c_schema["test_genes"])
    if train_genes & val_genes:
        add_error(errors, f"C train and validation gene lists overlap: {sorted(train_genes & val_genes)}")
    if train_genes & test_genes:
        add_error(errors, f"C train and test gene lists overlap: {sorted(train_genes & test_genes)}")
    if val_genes & test_genes:
        add_error(errors, f"C validation and test gene lists overlap: {sorted(val_genes & test_genes)}")
    gene_split = {}
    for gene in c_schema["train_genes"]:
        gene_split[gene] = "train"
    for gene in c_schema["validation_genes"]:
        gene_split[gene] = "validation"
    for gene in c_schema["test_genes"]:
        gene_split[gene] = "test"
    for message in gene_partition_errors(gene_split):
        add_error(errors, f"C_heldout_target_gene partition: {message}")
    c_rows = [
        (row.condition, row.C_heldout_target_gene)
        for row in split_obs.itertuples(index=False)
    ]
    for message in target_gene_split_errors(c_rows, gene_split):
        add_error(errors, f"C_heldout_target_gene: {message}")

    d_map = {int(k): v for k, v in proposal["protocols"]["D_heldout_barcode_group"]["group_to_split"].items()}
    bad_d = split_obs.loc[
        split_obs["barcode_group"].map(d_map) != split_obs["D_heldout_barcode_group"]
    ]
    if len(bad_d):
        add_error(
            errors,
            f"D_heldout_barcode_group has {len(bad_d)} rows that do not follow group_to_split",
        )

    if evidence_path.is_file():
        evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
        if evidence.get("dataset_sha256") != expected_sha:
            add_error(
                errors,
                "provenance_evidence.json dataset_sha256 does not match schema identity.sha256",
            )
        if evidence.get("target_genes") != expected_targets:
            add_error(errors, "provenance_evidence.json target_genes != schema targets.names")
        json_splits = evidence.get("splits", {})
        for protocol, payload in proposal["protocols"].items():
            for split_name, expected in payload["sizes"].items():
                actual = json_splits.get(protocol, {}).get(split_name, {})
                if actual.get("cells") != expected["cells"]:
                    add_error(
                        errors,
                        f"JSON {protocol}.{split_name}.cells={actual.get('cells')} "
                        f"!= schema {expected['cells']}",
                    )
    else:
        add_error(errors, f"missing provenance evidence: {evidence_path}")

    guide_path = ROOT / schema["source_provenance"]["verified_guide_sequences"]
    if guide_path.is_file():
        guides = pd.read_csv(guide_path)
        if "bulk_fitness" in guides.columns:
            add_error(errors, "verified guide table must not contain bulk_fitness")
        if set(guides["stored_target_gene"]) != set(expected_targets):
            add_error(errors, "verified guide table target names != schema 105-name list")
        unresolved = guides.loc[guides["unresolved_status"].astype(str).str.startswith("unresolved")]
        expected_unresolved = set(schema["source_provenance"]["unresolved_guide_targets"])
        if set(unresolved["stored_target_gene"]) != expected_unresolved:
            add_error(
                errors,
                f"unresolved guide targets {sorted(unresolved['stored_target_gene'])} "
                f"!= schema {sorted(expected_unresolved)}",
            )
    else:
        add_error(errors, f"missing verified guide-sequence table: {guide_path}")

    return errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--schema",
        type=Path,
        default=DEFAULT_SCHEMA,
        help="Path to dataset_schema.yaml",
    )
    args = parser.parse_args()
    errors = validate(args.schema)
    if errors:
        print(f"SCHEMA VALIDATION FAILED ({len(errors)} error(s))", flush=True)
        for message in errors:
            print(f"ERROR: {message}", flush=True)
        return 1
    print("SCHEMA VALIDATION PASSED", flush=True)
    print(f"schema={args.schema}", flush=True)
    print("training_allowed=false final_protocol_approved=false", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
