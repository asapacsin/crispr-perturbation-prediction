## Audited biological summary

Numbers below are copied from [MODELING_DIAGNOSTICS.md](../MODELING_DIAGNOSTICS.md) and `reports/modeling/diagnostics/`. They did not enter fit or alpha selection.

**Plain-language example.** On stored log1p `X`, CEBPA alone raises LST1 by 2.46883, CEBPB alone by 1.14186, and the pair reaches 2.22471 rather than the additive 3.61069. The LST1 residual is -1.38597 (conditional SE 0.165184). Root recomputed that same LST1 contrast on original stored log1p `X` ([root_verification.json](root_verification.json)). This is one-experiment cell-sampling uncertainty, not a biological p-value. The pair has **n=75** cells. Construct-orientation alternative additive RMSE for this pair is **0.09046–0.12455**, so a residual remains but its magnitude is orientation-sensitive. Do not read this as proven synergy.

Top 3 |residual| genes in `CEBPA+CEBPB`:

| gene | A | B | AB | A+B | residual | SE |
| --- | --- | --- | --- | --- | --- | --- |
| LST1 | 2.46883 | 1.14186 | 2.22471 | 3.61069 | -1.38597 | 0.165184 |
| AIF1 | 1.44996 | 0.998837 | 1.16211 | 2.44879 | -1.28668 | 0.137779 |
| CFD | 2.10394 | 0.959009 | 1.93222 | 3.06295 | -1.13073 | 0.195137 |

**Own-target check.** 103 of 105 measured targets have a positive own-locus log1p effect (2 negative). That is a normalized consistency check, not molecular validation of guides or protein.

**A more additive contrast.** `CEBPA+ZC3HAV1` has residual/response RMS **0.275** (observed RMS 0.165808, residual RMS 0.0456002). That is closer to additivity than CEBPA+CEBPB and is **not** perfect additivity.

**Modules, not new GO.** Named GO enrichment is unavailable (graph has source/target/importance only). Data-driven SVD PC1 separates LST1/AIF1/CSF3R from HBZ/HBG1/HBG2. Norman et al. 2019 already interpret related programs as granulocyte- versus erythroid-like; that is a **supported existing interpretation**, not a new enrichment test.

**Shared shift.** Median residual squared-norm fraction along the mean training single-delta is **0.124**. Alignment does not exclude other confounding.

**Held-out prediction (already frozen).** Ridge α=0.1 test macro MSE **0.001599794** vs additive **0.002372021** (32.56% lower; 15/26 pairs). The 26 test pairs are consumed. Next models need nested resampling or new external holdout.

Full diagnostic report: [MODELING_DIAGNOSTICS.md](../MODELING_DIAGNOSTICS.md). Figures: [pair ranking](describe/figures/pair_error_distribution_ranking.png), [zero-centered interaction heatmap](describe/figures/strongest_interaction_heatmap.png).

