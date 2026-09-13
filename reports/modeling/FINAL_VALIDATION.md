# Final validation

## Status: complete

Root executed the remaining packaging command. No `--stage fit`. Matplotlib **3.11.2**. Original H5 unchanged; `counts` unused.

| Check | Result |
| --- | --- |
| Unittests (close-out) | **71 PASS** |
| Schema validation | PASS |
| Dataset SHA-256 | `23ffb0fac6a847ff927cf7509d80d85052bfefbfb97610786a2dafaaefa0b6a0` |
| `--plots-only` | PASS |
| `finalize_modeling_report.py` | PASS |
| Pair-ranking PNG | Root inspected: labels clear, no overlap |
| Ridge α | 0.1 (validation-selected; not retuned) |
| `root_verification.json` | original stored log1p `X` example matches export; `all_metric_checks_pass: true` |
| `diagnostics/metric_verification.json` | `all_checked_metrics_match: true`; train ridge shape [184, 5045] |

Executed command:

```powershell
.\.venv-modeling\Scripts\python.exe -m unittest discover -s tests -v; .\.venv-modeling\Scripts\python.exe scripts\validate_dataset_schema.py; .\.venv-modeling\Scripts\python.exe scripts\modeling_analysis.py --plots-only; .\.venv-modeling\Scripts\python.exe scripts\finalize_modeling_report.py
```

## History: packaging-pass worker shell (blocked)

A worker attempted the same integrity loop and received `Rejected:` (Cursor shell permission denied). No bypass or force was used. No hash or unittest count was produced in that pass. Root later ran the command above successfully.

Exact blocked command (history only):

```powershell
.\.venv-modeling\Scripts\python.exe -m unittest discover -s tests -v; .\.venv-modeling\Scripts\python.exe scripts\validate_dataset_schema.py; .\.venv-modeling\Scripts\python.exe scripts\modeling_analysis.py --plots-only; .\.venv-modeling\Scripts\python.exe scripts\finalize_modeling_report.py; .\.venv-modeling\Scripts\python.exe -c "import hashlib,sys; from pathlib import Path; p=Path('data/norman/perturb_processed.h5ad'); h=hashlib.sha256(); f=p.open('rb');
[h.update(c) for c in iter(lambda: f.read(1<<20), b'')]; f.close(); print(h.hexdigest())"
```
