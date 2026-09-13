# Repository Guidelines

This project audits Norman / GEARS-processed Perturb-seq data and now includes an **authorized Package-B linear baseline pipeline**. Do not modify `data/norman/*`. Do not rerun `scripts/download_norman.py` casually: saved download headers can overwrite the dataset.

Current modeling authorization is `config/modeling_protocol.yaml` (lane-balanced effects + linear baselines; ridge α already selected at 0.1). `config/dataset_schema.yaml` `final_protocol_approved` / `training_allowed` remain false as **historical proposal flags**, not a denial of the completed baseline run. Do not refit or retune on the consumed 26 test pairs.

## Layout

- `scripts/`: inspection, validation, provenance, modeling, diagnostics, and `finalize_modeling_report.py`.
- `config/dataset_schema.yaml`: read-only H5AD identity plus proposal-only split fields.
- `config/modeling_protocol.yaml`: authorized modeling protocol.
- `data/norman/`: archive and extracted `perturb_processed.h5ad` (immutable).
- `reports/`: audits, provenance, `MODELING_REPORT.md`, `MODELING_DIAGNOSTICS.md`, `RESULTS.md`, and `reports/modeling/` artifacts.
- `tests/`: `unittest` modules (`test_dataset_schema.py`, `test_modeling.py`, `test_modeling_diagnostics.py`, `test_finalize_report.py`).

## Commands

Prefer the isolated venv created by root (matplotlib 3.11.2; no raw or base-env changes):

```powershell
.\.venv-modeling\Scripts\python.exe scripts\validate_dataset_schema.py
.\.venv-modeling\Scripts\python.exe -m unittest discover -s tests -v
.\.venv-modeling\Scripts\python.exe scripts\modeling_analysis.py --plots-only
.\.venv-modeling\Scripts\python.exe scripts\finalize_modeling_report.py
```

Historical audit/provenance scripts (do not rerun downloads):

```powershell
.\.venv-modeling\Scripts\python.exe scripts\validate_dataset_schema.py
.\.venv-modeling\Scripts\python.exe scripts\write_verified_guide_sequences.py
```

`validate_dataset_schema.py` hashes the H5AD and decodes obs/var categoricals without loading dense `counts`. Tests use tiny HDF5 fixtures and do not train models. `--stage fit` must not be rerun for packaging.

## Style & Tests

Use four-space indentation, `snake_case`, `pathlib`, and UTF-8. Name tests `tests/test_*.py`; cover identity, leakage, modeling math, and diagnostics.

## Commits & Pull Requests

The project root has no Git history to establish a convention. Use imperative subjects, such as `Validate perturbation schema`. PRs should state the change, evidence, validation commands, and unresolved limitations.

## Integrity

Preserve stored identifiers and barcode suffixes. Distinguish stored labels from resolved provenance. Existing `uns` DE rankings are prohibited for fitting. Keep large datasets, temporary downloads, and credentials out of commits.
