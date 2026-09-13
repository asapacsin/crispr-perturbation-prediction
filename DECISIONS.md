# Decisions

Scope decisions for the authorized modeling work. Historical schema and provenance files were not rewritten.

## Binding

1. **X is stored log1p.** Load sparse `X` via h5py CSR (backed AnnData only as a fallback). Never load the dense `counts` layer. Never read `uns`. Feature order is the stored `gene_id` order (5,045 genes).
2. **Condition unit.** Canonical unordered target sets group cells. Original orientation labels are preserved. Pooling `GENE+ctrl` with `ctrl+GENE` is a construct-identity caveat, not experimental equivalence.
3. **Lane-balanced means.** A condition mean is the average of the eight within-lane means. Single effect in a lane is `mean(single,lane)-mean(ctrl,lane)`. Pair interaction in a lane is `mean(pair,lane)-mean(A,lane)-mean(B,lane)+mean(ctrl,lane)`, then lanes are averaged. This is additivity on log1p space, not mechanistic synergy.
4. **Uncertainty.** Lane-mean variance is sample variance / n. Contrast variances are summed; the SE of the eight-lane average is `sqrt(sum / 64)`. Control enters a contrast once. If any component lane has n<2, SE is NA and flagged. No biological p-values or FDR. Significance field is `not_estimable_independent_bioreplicates_absent`. Eight lanes are not eight bioreplicates. Intervals are **conditional on this one experiment**.
5. **Order.** Predictive fit and evaluation happen first. All-data describe cannot change the protocol or the chosen ridge alpha.
6. **Split.** Existing manifest column `B_heldout_combination`, seed 20260911: 105 singles + 79 pairs train, 26 / 26 pairs val/test. Train-only transforms, centroids, and `mu0`. Val/test controls are diagnostic only.
7. **Baselines.** Control-only; train condition-balanced perturbed mean; matching-single mean; additive `mu_A+mu_B-mu0`; ridge multi-hot of 105 genes → 5,045-gene deltas, unpenalized intercept, centered design/response. Alphas `[0.1, 1, 10, 100, 1000]`; one alpha by validation macro MSE; no train+val refit; test once.
8. **Metrics.** Primary: condition-macro MSE. Secondary: MAE. **Pearson is NA if either vector is constant. Cosine is NA if either vector has zero norm.** Skills vs perturbed mean and additive; NA if baseline loss is 0. Paired bootstrap of test-condition MSE vs additive (seed 20260913, 2000) is conditional on the frozen model and train reference, not a biological CI.
9. **Schema flags.** `config/dataset_schema.yaml` stays unapproved. Current authorization is `config/modeling_protocol.yaml`.
10. **Test consumed.** The 26 test pairs have been evaluated and inspected. New models need a **new honest evaluation** (nested condition-level resampling as a retrospective benchmark, or new external holdout). Do not retune on these 26 pairs.
11. **No deep model from this scorecard.** Named GO enrichment is skipped. Diagnostics are descriptive.

## Explicitly not decided by this work

Whether orientation labels are interchangeable; whether residuals are mechanistic synergy; donor-level uncertainty; GEO UMI reconstruction; `SLC38A2` sequence; RHOXF2 family identity; Phase 2 GEARS or any distributional/deep package.
