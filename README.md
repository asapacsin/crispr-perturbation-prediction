# CRISPR perturbation prediction

Audit of Norman / GEARS-processed Perturb-seq plus an authorized Package-B linear baseline. The target is lane-balanced log1p effects (91,205 cells × 5,045 genes; K562 CRISPRa). A ridge multi-hot main-effect map is compared with additive prediction on held-out combinations. No deep model was trained.

## Results

Ridge α **0.1** was selected on validation and evaluated once. Test macro MSE: ridge **0.001599794** vs additive **0.002372021** (32.56% lower; 15/26 pairs). Own-target sign is positive for 103/105 singles. The 26 B-protocol test pairs are consumed; do not retune on them.

Full numbers and limitations: [RESULTS.md](RESULTS.md). Methods, tables, and figures: [reports/MODELING_REPORT.md](reports/MODELING_REPORT.md).

## Data (not in this repository)

The processed H5AD is excluded. Download the Dataverse archive that contains `norman/perturb_processed.h5ad`:

- Archive: https://dataverse.harvard.edu/api/access/datafile/6154020
- Dataset: https://doi.org/10.7910/DVN/Q2ZV3E

Place the extracted file at `data/norman/perturb_processed.h5ad`. Expected SHA-256: `23ffb0fac6a847ff927cf7509d80d85052bfefbfb97610786a2dafaaefa0b6a0`. Do not rerun `scripts/download_norman.py` casually (saved headers can overwrite the file). Do not load dense `counts` for schema checks.

## Reproduce (no refit)

Create a local environment (Python 3.14 was used for the recorded run):

```powershell
python -m venv .venv-modeling
.\.venv-modeling\Scripts\python.exe -m pip install -r requirements.txt
```

The commands below validate and render the saved results without refitting.

```powershell
.\.venv-modeling\Scripts\python.exe scripts\validate_dataset_schema.py
.\.venv-modeling\Scripts\python.exe -m unittest discover -s tests -v
.\.venv-modeling\Scripts\python.exe scripts\modeling_analysis.py --plots-only
.\.venv-modeling\Scripts\python.exe scripts\finalize_modeling_report.py
```

Protocol: `config/modeling_protocol.yaml`. Identity and historical proposal flags: `config/dataset_schema.yaml`.
