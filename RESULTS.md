# RESULTS

Compact final package. **Executed** = observed by root. Historical worker shell block is history only.

## 1. Objective

Estimate log1p single/pair effects on Norman Perturb-seq and test whether a linear main-effect map beats additive prediction on held-out combinations (Package B). Priority: interpretability, not a deep model.

## 2. Methods actually used

Stored CSR `X` (log1p; no dense `counts`, no `uns`). Lane-balanced means (8 technical GEM lanes). Single effect and pair interaction as predeclared. Split: existing `B_heldout_combination` seed 20260911 (105 singles + 79 pairs train; 26/26 val/test). `mu0` = train controls only. Baselines: control-only, train condition-balanced perturbed mean, matching-single, additive, ridge multi-hot 105→5045 with unpenalized intercept; α∈{0.1,1,10,100,1000} selected on val macro MSE, frozen, test once.

## 3. Implementation decisions

Authorization in `config/modeling_protocol.yaml`; schema approval flags left false. Pearson NA if constant; cosine NA if zero-norm. Diagnostics after fit only. Isolated `.venv-modeling` (matplotlib 3.11.2) created by root — not a global no-install rule. Ranking plot uses `constrained_layout`. Audited prose is generated from diagnostic tables by `finalize_modeling_report.py` with destination-relative links.

## 4. Validation (executed)

- Fit-era: 69 tests PASS and schema PASS (root, before packaging).
- Close-out (root): **71 unittests PASS** (includes `test_finalize_report.py`), schema PASS, H5 SHA-256 `23ffb0fac6a847ff927cf7509d80d85052bfefbfb97610786a2dafaaefa0b6a0`, `--plots-only` PASS, finalizer PASS. No fit rerun. Matplotlib 3.11.2.
- `metric_verification.json`: `all_checked_metrics_match: true`.
- `root_verification.json`: original stored log1p `X` CEBPA+CEBPB / LST1 matches export; independent test MSEs match.

## 5. Key numbers (executed; no refit)

| Quantity | Value |
| --- | --- |
| Chosen α (val MSE) | 0.1 (0.00137707) |
| Test MSE ridge / additive | 0.001599794 / 0.002372021 (32.56% lower; 15/26 pairs) |
| Bootstrap MSE(ridge−additive) | −0.000772; 2.5/97.5: −0.001449 / −0.000224 (conditional) |
| Own-target sign | 103 / 105 positive |
| CEBPA+CEBPB n | 75 cells (test pair) |
| LST1 A / B / AB / residual / SE | 2.468826 / 1.141860 / 2.224714 / −1.385972 / 0.165184 |
| Construct alt. RMSE that pair | 0.09046–0.12455 |
| CEBPA+ZC3HAV1 resid/response | 0.275 (not perfect additivity) |
| Median residual‖proj_c‖²/‖r‖² | 0.124 |

## 6. Unexpected / resolved engineering

Initial schema reader failed on legacy `__categories` datasets; fixed. NumPy scalar-test / `json_safe` bool-before-int bugs fixed. Ranking-figure y-labels overlapped the left histogram; `constrained_layout` helper added; root regenerated PNG and confirmed labels clear, no overlap.

## 7. Limitations

One pooled K562 CRISPRa experiment; 8 lanes ≠ bioreplicates. Full-cohort HVG. Additivity is log1p arithmetic. CI is conditional cell-sampling. 26 test pairs consumed. No biological p-values. No deep-model claim.

## 8. Unresolved scientific questions

Orientation interchangeability; whether residuals are synergy vs construct/composition; donor-level uncertainty; GEO UMI rebuild; SLC38A2 / RHOXF2 identity.

## 9. Important paths

`reports/MODELING_REPORT.md`, `reports/MODELING_DIAGNOSTICS.md`, `reports/modeling/AUDITED_SUMMARY.md`, `reports/modeling/root_verification.json`, `reports/modeling/diagnostics/metric_verification.json`, `reports/modeling/fit/summary_metrics.json`, `reports/modeling/describe/figures/pair_error_distribution_ranking.png`, `reports/modeling/describe/figures/strongest_interaction_heatmap.png`, `reports/modeling/FINAL_VALIDATION.md`.

## 10. Packaging history

**Worker (files only; shell `Rejected:`):** report splice, helpers, `test_finalize_report.py`. Ranking PNG then still overlapped. Integrity commands were not run from that worker.

**Root (close-out command executed):** 71 tests, schema+SHA, `--plots-only`, finalizer. Ranking labels inspected clear. Original H5 unchanged; `counts` unused.
