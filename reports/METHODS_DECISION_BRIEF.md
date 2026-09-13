# Methods decision brief — Norman / GEARS-processed Perturb-seq

**Decision support only. No task, split, metric, method, visualization, or budget is approved.**
Revised 2026-09-12 after reading three retrieved papers. No model trained; H5AD unchanged.

## Local facts (binding)

| Fact | Evidence | Protocol implication |
| --- | --- | --- |
| 91,205 × 5,045; `X` = per-cell scaled natural log1p | `DATA_AUDIT.md` | Estimand on stored `X`, not log of aggregated counts |
| 5,045 genes retained on the **full** object | same | Not fully inductive; reprocessing needed for train-only HVGs |
| 105 singles + 131 pairs + control; 7,353 `ctrl` | same | Canonical target sets as the condition unit |
| 8 GEM lanes of one pooled K562 CRISPRa experiment | GEO + `barcode_group_summary.csv` | Technical partitions, not 8 biological replicates |
| Stored `cell_type` = `A549` | GEARS issue 95; preprocessing snapshot | Dummy field; do not model A549 |
| No matched before/after cells | design | No individual causal effect |
| `uns` DE on the full object | `DATA_AUDIT.md` | Do not use `full.uns` in fitting/tuning |

Proposed splits (`provenance_evidence.json`, seed 20260911) remain **unapproved**:

| Plan | Train (ctrl) | Val (ctrl) | Test (ctrl) | Hold-out |
| --- | ---: | ---: | ---: | --- |
| A | 73,219 (5,883) | 8,993 (735) | 8,993 (735) | Cells only; all conditions in every split |
| B | 75,588 (5,891) | 7,639 (731) | 7,978 (731) | 105 singles + 79 pairs train; 26 / 26 pairs val/test |
| C | 59,007 (5,891) | 12,676 (731) | 18,476 (731) | 73 / 16 / 16 genes; 4 bridge pairs unused |
| D | 67,801 (5,470) | 11,667 (934) | 11,737 (949) | Lanes 1–6 / 7 / 8 |

C bridges (1,046 cells): `FOXA1+FOXF1`, `MAPK1+TGFBR2`, `BAK1+BCL2L11`, `BAK1+TMSB4X`.

---

## 1. Primary-source evidence

No ranking. Derived prose per paper is kept short. GEARS / CPA / biolord / PerturBench publisher pages remain **unretrieved**.

| Source | Access | Establishes | Does **not** establish |
| --- | --- | --- | --- |
| [GEARS](https://www.nature.com/articles/s41587-023-01905-6) Roohani, Huang, Leskovec, *Nat. Biotechnol.* 2023/24 | Paper **unretrieved**. Local repo v0.1.2 read | Graph method (GO GNN + co-expression GNN) for novel single/multi-gene outcomes. Norman loader uses Dataverse `6154020`. README: “When trained on single-gene perturbation data alone, GEARS cannot reliably predict outcomes for combinatorial perturbations.” | Wins on these splits; combo prediction from singles only |
| [Ahlmann-Eltze, Huber & Anders](https://www.nature.com/articles/s41592-025-02772-6) *Nat. Methods* **22**, 1657–1661 (2025); published 4 Aug 2025. Local `methods_linear_benchmark.txt` | **Retrieved** | On Norman doubles, no tested deep/foundation model beat **no-change** or **additive** (L2 on 1,000 highly expressed genes; five splits). On unseen singles (Adamson/Replogle), none consistently beat a **mean** or **bilinear/ridge** baseline. Interactions vs additive were no better than no-change. | That no future deep model can win; a score on *this* 91,205 × 5,045 file |
| [Viñas Torné et al.](https://www.nature.com/articles/s41587-025-02777-8) *Nat. Biotechnol.* **44**, 1050–1059 (2026); published 25 Aug 2025. Local `methods_systema.txt` | **Retrieved** | **Systematic variation** (shared perturbed-vs-control shift) inflates Pearson-Δ. Their **perturbed mean** is the **cell-weighted** average of all train perturbed cells; **matching mean** averages seen-partner centroids (unseen partner → perturbed mean). Norman is high-systematic-variation. **MSE/RMSE is reference-insensitive**: subtracting the same vector from truth and prediction does not change it. Systema recentres on a **perturbed-centroid** to score specific effects / centroid accuracy. | A universal winner; that cell-weighted mean equals a condition-balanced mean |
| [Wong, Hill & Moccia](https://pmc.ncbi.nlm.nih.gov/articles/PMC12202205/) *Bioinformatics* **41**, btaf317 (2025); PMC12202205. Local `methods_simple_controls.txt` | **Retrieved** | Same Norman Dataverse `6154020`. A train **perturbed-cell mean** was competitive with GEARS/scGPT. Their stronger “CRISPR-informed mean” **zeros** the target under CRISPRi and **doubles** it under CRISPRa — “an overly simplistic approximation.” On Norman CRISPRa, that informed mean did **not** dominate DE-delta vs GEARS. Perturb-seq is unpaired. | That CRISPRi zeroing is appropriate here; that mean captures gene–gene interaction |
| [CPA](https://doi.org/10.15252/msb.202211517); [biolord](https://doi.org/10.1038/s41587-023-02079-x); [PerturBench](https://www.nature.com/articles/s41592-025-02980-0) | **Unretrieved** | Assigned roles only: optional distribution/disentanglement comparators; recent benchmark hygiene. No leader claimed | Any ranking on this file |
| Norman et al. 2019 + GEO | Local PMC + `GSE133344` / `GSM3906020` | One pooled K562 CRISPRa Perturb-seq; 8 Chromium lanes from “the same experiment.” Source GI model: linear combination of single deltas. “K562 cells stably expressing the SunTag CRISPRa system.” | Paired cells; that 91,205 is the complete original call set; that GEO’s `dCas9-KRAB` line overrides the paper |

A549 is a dummy (`GEARS_issue95_comments.json`; preprocessing sets `cell_type = 'A549'`).

**Leakage vs legitimate prior.** Using **external** GO membership, sequences, or other descriptors of **known gene identities** (including test gene *names* as graph nodes) is a **transductive / external prior**, not outcome leakage. Leakage is using **held-out expression or outcomes**: co-expression edges from test cells, `full.uns` DE, test-derived HVGs, or any graph rebuilt from val/test `X`.

---

## 2. Candidate methods

Costs are engineering estimates. **One-hot ridge cannot extrapolate to completely unseen target identities** without informative external descriptors.

**Condition-balanced perturbed mean** (this proposal): unweighted mean of **train condition centroids** (optionally average singles and pairs separately). It is **always defined** whenever any train perturbed conditions exist (true for A–D). It is **not** Systema’s cell-weighted perturbed mean, which weights large conditions more.

| Method | Predicts | Assumptions | B | C | Risk | Cost |
| --- | --- | --- | --- | --- | --- | --- |
| Train-only no-change | Train control mean | Effects small vs control | Required null | Required null | Freeze train controls | CPU, min |
| Condition-balanced perturbed mean | Equal-weight mean of train condition centroids | Shared shift across conditions | Defined | Defined (common shift only) | Do not use test to form it | CPU, min |
| Matching single mean | Train mean of *i*, *j*, or average | One partner dominates | Defined (both singles in train) | Only for the seen partner; else fall back to balanced mean (pre-specify) | Freeze partner rule | CPU, min |
| Additive deltas | \(\mu_{\mathrm{ctrl}}^{\mathrm{train}}+\delta_i+\delta_j\) | Log1p effects add | Primary biological baseline | Undefined for unseen genes | Freeze orientation collapse | CPU, min |
| Regularized low-rank / ridge | Linear map from train-visible encodings | Residual ≈ low-rank | Linear comparator | Needs external descriptors, not one-hot IDs | Medium if test DE used as features | CPU, min–1 h |
| GEARS | Graph-informed transcriptome | External GO ± **train-only** co-expression | Optional after additive/ridge; needs some train pairs | Eligible via GO prior, not a winner; genes off the GO graph are dropped | High: PyG, split mismatch, their `_de` slice | GPU hours |
| CPA (only if Package 3) | Cell-level distribution | Estimand is a distribution | Not warranted for a mean task | Same | High | GPU hours–day |

biolord is not a second default distribution model.

**Phase-1 proposal (not approval):** no-change, condition-balanced perturbed mean, matching-single (B), additive (B), optional ridge. GEARS only if Phase 2 is funded. CPA only if Package 3 is chosen.

---

## 3. Estimand (proposed, not approved)

For canonical target set \(c\): condition mean of stored log1p, and/or its delta to the **train-only** control mean.

\[
\mu_c=\mathbb{E}[X_i\mid i\in c],\qquad \delta_c=\mu_c-\mu_{\mathrm{ctrl}}^{\mathrm{train}}
\]

Not \(\log(1+\sum\mathrm{count})\), not a paired individual effect, not count-pseudobulk DE (mean-log ≠ aggregate-then-test).

**Same-reference MSE identity.** For any shared vector \(r\) (control mean or common shift), \(\mathrm{MSE}(\mu-r,\hat\mu-r)=\mathrm{MSE}(\mu,\hat\mu)\). Do **not** claim that “delta MSE” differs numerically from absolute MSE when the same \(r\) is subtracted from both sides. MAE is likewise unchanged by a shared additive shift.

**Readout (choose one before tuning).**

| Option | Meaning | Cost of the choice |
| --- | --- | --- |
| Full 5,045 | All stored genes | Full-object filter |
| Target-excluded 4,940 | Drops **all 105** assigned targets, not only the one or two genes in the queried condition | Removes the overexpressed loci from the score, but also drops **other panel TFs/effectors** that are often the real downstream program in this CRISPRa design (e.g. scoring KLF1/CEBP/hemoglobin-related targets when they are not the intervention). Still not an inductive HVG set |
| Train-only features | Re-filter from a broader count matrix on train cells | Requires upstream reprocessing this H5AD does not support |

Prepared data are **not** fully leakage-free (5,045 filter; `uns` DE; any expression graph that includes held-out `X`).

---

## 4. Metrics (choose before seeing test)

**One primary (proposal, not approval): condition-macro MSE** of \(\hat\mu_c\) vs \(\mu_c\) (equivalently of deltas vs the same train control). Equal weight per test condition.

- **Secondary:** condition-macro **MAE**.
- **Skill (report, not a second primary):** \(1-L(\hat\mu)/L(\text{baseline})\) vs **(i)** condition-balanced perturbed mean and **(ii)** additive on B. No-change may be reported but is not the skill pair that justifies extra methods once a common shift exists.
- **Pearson of deltas:** secondary only. If either vector is constant, record `NA`; do not coerce to 0. Subtracting a **shared** common shift **does not change MSE/MAE**; it only changes correlation/direction. To ask “beyond systematic shift,” compare **skill to the perturbed-mean baseline**, or inspect residual **direction** after a **train-only** common shift — do not recentre on test.

Do not fit/tune on `full.uns` DE lists. Test-derived DE panels are post-freeze secondary only (winner-selection and noise caveats).

Bootstrap **conditions**; optional within-condition cell resampling. Eight lanes ≠ eight biological replicates.

---

## 5. Three task packages

**B — held-out combinations (proposed, not approved).** Predict mean log1p of 26/26 unseen pairs after 105 singles + 79 pairs in train. This identifies **predictive combination generalization** (can a map from seen singles/pairs forecast a new pair’s mean?). The split alone does **not** identify a causal interaction or a paired treatment effect. Methods: no-change, balanced perturbed mean, matching single, additive, ridge, optional GEARS. “Beat additive” is defined only here.

**C — held-out target genes.** 73 / 16 / 16 genes; exclude the four bridges. Identifies transfer to new **identities**. One-hot ridge and additive are not applicable. External GO/descriptors are legitimate priors; expression-from-test graphs are not.

**3 — cell distributions (only if chosen).** Same splits, distributional estimand (e.g. sliced Wasserstein / MMD / gene marginals). **Required baselines:** (1) resample train **controls** as a null distribution; (2) empirical distribution of the **matching single** (B) or of the **condition-balanced train perturbed cells** — not only a degenerate point mass at a mean. CPA is then the one deep distribution model. Still no paired counterfactuals. Choose this *instead of* the mean task.

**Attached, not packages:** D = lane transport; A = interpolation sanity.

---

## 6. Visualization packages (do not generate)

Freeze with the task; no test-outcome cherry-picking.

| Package | Answers | Hides | Status |
| --- | --- | --- | --- |
| **V1 recommended** | Scorecard (primary MSE, MAE, skill vs perturbed-mean and additive) + paired per-condition MSE comparisons against each baseline + one residual heatmap with genes/conditions listed **before** test look | Heterogeneity; unlisted genes | Standard scorecard/error plots and heatmaps; exploratory until frozen, final reporting after results and a frozen list |
| **V2** | Absolute-expression UMAP / Pearson of \(\mu\) | Shared shifts; makes any mean look strong | Common in single-cell papers; exploratory only |
| **V3** | DE-gene / GI gallery | Full-transcriptome skill; leaks if genes from `full.uns` | Common in perturbation papers; post-freeze secondary |

---

## 7. Budget and decision boundary

Do not train a final system before the user approves **task/split**, **the single primary metric + gene set**, and **method budget**.

| Phase | Work | Stop |
| --- | --- | --- |
| 0 | User chooses B / C / 3, optional D/A, estimand, V1 | No fitting |
| 1 (default proposal, not approval) | No-deep set above | If perturbed-mean (and additive on B) already meet the bar, stop |
| 2 | GEARS only if funded or a declared additive gap on B remains | Frozen primary MSE; no `uns` DE in tuning |
| 3 | CPA only if Package 3 approved, vs control-resample and empirical-single distributions | No extra distribution model by default |

**Explicit default proposal, not approval:** Package **B** + stored-log1p condition mean + **condition-macro MSE** (MAE secondary; skill vs balanced perturbed-mean and additive) + V1 later + Phase-1 methods only.
