# Project state

Last updated 2026-09-13. Modeling analysis and packaging are **complete**. The 26 B-protocol test pairs are **consumed**.

## What this repository is

A read-only audit of the local Norman / GEARS-processed Perturb-seq H5AD plus an authorized Package-B linear baseline and interpretable effect analysis. `data/norman/*` is immutable. There is no deep-learning pipeline.

## Authorization

`config/modeling_protocol.yaml` authorizes lane-balanced effects and Package B linear baselines. Historical `config/dataset_schema.yaml` flags remain `final_protocol_approved: false` and `training_allowed: false` (not edited). Those flags are historical proposal records, not the current modeling authorization.

## Execution status

Root ran prepare → fit → describe and diagnostics. Ridge α **0.1** was selected on validation. Test macro MSE **0.001599794** vs additive **0.002372021**. Diagnostics completed (103/105 own-target positive; CEBPA+CEBPB LST1 residual verified on original stored log1p `X`). Original H5 unchanged; `counts` unused.

Close-out (root): **71 unittests PASS**, schema PASS (H5 SHA-256 `23ffb0fac6a847ff927cf7509d80d85052bfefbfb97610786a2dafaaefa0b6a0`), `--plots-only` PASS, finalizer PASS. Regenerated pair-ranking labels inspected clear (no overlap). Matplotlib 3.11.2. No fit rerun.

A packaging-pass worker shell was `Rejected:` (permission denied). That block is history; root later ran the remaining command. Isolate-venv execution uses `.venv-modeling` (matplotlib 3.11.2). Root created that venv; the base environment and stored data were not changed. A prior worker-only “no installs” note was **not** a global user policy.

## Commands

```powershell
.\.venv-modeling\Scripts\python.exe -m unittest discover -s tests -v
.\.venv-modeling\Scripts\python.exe scripts\validate_dataset_schema.py
.\.venv-modeling\Scripts\python.exe scripts\modeling_analysis.py --plots-only
.\.venv-modeling\Scripts\python.exe scripts\finalize_modeling_report.py
```

Do **not** rerun `--stage fit`. Do not reuse the 26 test pairs for tuning.

## Unchanged on purpose

`data/norman/*`, `config/dataset_schema.yaml`, historical provenance/audit reports, proposed split manifest.
