# Data provenance — Norman 2019 processed H5AD

This report joins the local processed file to primary sources. It does **not** replace [DATA_AUDIT.md](DATA_AUDIT.md), which remains the historical pre-provenance audit and must stay unchanged.

**Status:** provenance documented; schema validation is read-only. **No model was trained.** No proposed task, split, metric, model, or figure is approved. The user reserves those final choices. `final_protocol_approved: false` and `training_allowed: false` in [../config/dataset_schema.yaml](../config/dataset_schema.yaml).

**Not done here:** full GEO count-matrix reconstruction; new DE; plots; training; deletion or mutation of `data/norman/*`.

Evidence grades used below:

| Grade | Meaning |
| --- | --- |
| **Verified** | Directly checked against the local file and/or a downloaded primary source. |
| **Inferred** | Follows from verified facts plus an explicit rule; could still be wrong if the rule is. |
| **Unknown** | Not established by the files in this project. |

---

## 1. Claim-to-evidence table

| Claim | Grade | Primary source | Local evidence |
| --- | --- | --- | --- |
| Local H5AD SHA-256 `23ffb0fac6a847ff927cf7509d80d85052bfefbfb97610786a2dafaaefa0b6a0`, shape **91,205 × 5,045** | Verified | [Dataverse file 6154020](https://dataverse.harvard.edu/api/access/datafile/6154020); dataset [doi:10.7910/DVN/Q2ZV3E](https://doi.org/10.7910/DVN/Q2ZV3E) | [audit_evidence.json](audit_evidence.json), [provenance_evidence.json](provenance_evidence.json), `data/norman/perturb_processed.h5ad` |
| Only stored extra layer is `counts`; `X` is the main matrix, not a layer; no `raw` | Verified | File structure | [DATA_AUDIT.md](DATA_AUDIT.md), HDF5 keys `X, layers, obs, uns, var` |
| `X = log1p(counts × per-cell scale)` for all 37,317,477 nonzero entries; max abs error 6.13e-07 | Verified | Numerical identity on the local matrices | [validation_evidence.json](validation_evidence.json) |
| Numeric `normalize_total` target is the scanpy default (median nonzero library size), but the **number** is unknown and was not recovered from this H5AD | Inferred / unknown | [scanpy 1.8.1 `normalize_total`](https://github.com/scverse/scanpy/blob/1.8.1/scanpy/preprocessing/_normalization.py); GEARS notebook cell 25 | [scanpy181_normalization.py](provenance_sources/scanpy181_normalization.py), [Norman19_preprocessing_source.txt](provenance_sources/Norman19_preprocessing_source.txt) |
| All 91,205 barcodes join [GSE133344 raw cell identities](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133344) 1:1 | Verified | [GSE133344_raw_cell_identities.csv.gz](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE133344&format=file&file=GSE133344%5Fraw%5Fcell%5Fidentities.csv.gz) | [verified_cell_metadata.csv](verified_cell_metadata.csv), [provenance_evidence.json](provenance_evidence.json) |
| Joined rows are all `good_coverage=true` and `number_of_cells=1` | Verified | Same raw table | [provenance_evidence.json](provenance_evidence.json) (`source_good_coverage_counts`, `source_number_of_cells_counts`) |
| Raw good singlets 93,658 minus 2,453 `NegCtrl1_NegCtrl0__NegCtrl1_NegCtrl0` equals the local barcode set | Verified | Raw identities + notebook cell 13 | [provenance_evidence.json](provenance_evidence.json) (`raw_good_singlet_minus_excluded_equals_local_barcodes: true`); [Norman19_preprocessing_source.txt](provenance_sources/Norman19_preprocessing_source.txt) |
| 2,362 local cells are `cellranger_called=false` | Verified | Raw identities | [provenance_evidence.json](provenance_evidence.json) |
| Filtered GEO identities miss those 2,362 cells and conflict on 5 guide labels; **use raw** | Verified | [GSE133344_filtered_cell_identities.csv.gz](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE133344&format=file&file=GSE133344%5Ffiltered%5Fcell%5Fidentities.csv.gz) | [filtered_source_conflicts.csv](filtered_source_conflicts.csv) |
| One pooled Perturb-seq experiment; eight 10x Chromium **technical lanes / gemgroups**, not eight donors or biological replicates | Verified | [GSM3906020 SOFT](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3906020): “8 independent replicate samples”; “lanes on the 10X Chromium machine” | [GSM3906020.soft.txt](provenance_sources/GSM3906020.soft.txt); suffix = `gemgroup` for all 91,205 cells |
| Cell line is **K562** (high confidence) | Verified | Paper: “K562 cells stably expressing the SunTag CRISPRa system”; GEO sample source and characteristics name K562 | [Norman2019.html](provenance_sources/Norman2019.html) ¶P9–P10; [GSM3906020.soft.txt](provenance_sources/GSM3906020.soft.txt); [GSE133344.soft.txt](provenance_sources/GSE133344.soft.txt) |
| Stored `cell_type` / `condition_name` prefix `A549` is a **dummy** written by the GEARS preprocessor | Verified | Notebook cell 36: `adata_new.obs['cell_type'] = 'A549'`; GEARS maintainer: “That's just a dummy variable added by the dataloader.” | [Norman19_preprocessing_source.txt](provenance_sources/Norman19_preprocessing_source.txt); [GEARS issue 95](https://github.com/snap-stanford/GEARS/issues/95) |
| Mechanism is **SunTag CRISPRa overexpression**, not CRISPRi/KRAB knockdown | Verified paper/series; KRAB line is a **strong inference** of metadata error, not an author-confirmed correction | Paper ¶P10; series overall design describes overexpression. GSM3906020 also says “K562 cells expressing dCas9-KRAB”. Do **not** take that sentence as mechanism. | [Norman2019.html](provenance_sources/Norman2019.html); [GSE133344.soft.txt](provenance_sources/GSE133344.soft.txt) `!Series_overall_design`; [GSM3906020.soft.txt](provenance_sources/GSM3906020.soft.txt) `!Sample_growth_protocol_ch1` |
| Paper identity: Science 2019, DOI [10.1126/science.aax4438](https://doi.org/10.1126/science.aax4438), PMID [31395745](https://pubmed.ncbi.nlm.nih.gov/31395745/), PMC [PMC6746554](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6746554/) | Verified | Journal / PubMed / PMC / GEO `!Series_pubmed_id = 31395745` | [GSE133344.soft.txt](provenance_sources/GSE133344.soft.txt); [Norman2019.html](provenance_sources/Norman2019.html) |
| Preprocessor selected 5,000 HVGs plus targets → 5,045 genes; `normalize_total` default then `log1p` on **24,665** genes before subset | Verified | Notebook cells 25–31; notebook outputs `93658 × 24665` then `91205 × 24665` then `91205 × 5045` | [Norman19_preprocessing.ipynb](provenance_sources/Norman19_preprocessing.ipynb); commit [f882118](https://github.com/yhr91/GEARS_misc/blob/f88211870dfa89c38a2eedbd69ca1abd28a25f3c/data/preprocessing/Norman19.ipynb) |
| Full-data HVG selection is a leakage caveat for prospective claims; HVGs may still be used in a **retrospective prepared-feature benchmark**. Stored `uns` DE lists are **prohibited for fitting** | Inferred | HVG and `rank_genes_groups_by_cov` run on all retained cells before any of the proposed splits | Notebook cells 25, 43; [validation_evidence.json](validation_evidence.json) `uns` |
| Raw source label `RHOXF2` vs stored `RHOXF2BB` accounts for **762** rows; the other **90,443** labels agree before that mapping | Verified | Raw `guide_identity` vs stored `condition` | [provenance_evidence.json](provenance_evidence.json); [guide_condition_mapping.csv](guide_condition_mapping.csv) |
| Author relabel: guide lands between RHOXF2 and RHOXF2B, stronger reported effect on B | Verified as an author comment, not as a new molecular assay | [GI_generate_populations.ipynb](https://github.com/thomasmaxwellnorman/Perturbseq_GI/blob/3b25109aeb9c0c2026bd70abd50304a0ad4e5395/GI_generate_populations.ipynb) cell 3 | [GI_generate_populations.ipynb](provenance_sources/GI_generate_populations.ipynb) |
| 105 stored target tokens; 152 single labels → 105 target sets; 131 doubles; cells 7,353 / 48,407 / 35,445 | Verified | Stored `condition` | [perturbation_evidence.json](perturbation_evidence.json), [DATA_AUDIT.md](DATA_AUDIT.md) |
| 289 source `guide_identity` values collapse to 284 conditions after suffix variants and control merge. Never count both `__` halves as four perturbations | Verified | Parser uses the second construct half only | [guide_condition_mapping.csv](guide_condition_mapping.csv) (`guide_variants_per_condition` 280/3/1) |
| Author guide tables differ (105 vs 104 rows). `SLC38A2` has no source sequence. `RHOXF2BB` maps through documented `RHOXF2` | Verified | [data_sharing table](https://github.com/thomasmaxwellnorman/Perturbseq_GI/blob/3b25109aeb9c0c2026bd70abd50304a0ad4e5395/data_sharing/20180821_guides_chosen_for_perturbseq.csv), [final_emaps table](https://github.com/thomasmaxwellnorman/Perturbseq_GI/blob/3b25109aeb9c0c2026bd70abd50304a0ad4e5395/final_emaps/20180821_guides_chosen_for_perturbseq.csv) | [verified_target_guide_sequences.csv](verified_target_guide_sequences.csv); preliminary [target_guide_sequences.csv](target_guide_sequences.csv) is incomplete |
| `bulk_fitness` is an outcome-derived external score and must not be a model input | Verified as a leakage rule | Column exists only in the data_sharing table | Do not copy it into training features |

---

## 2. Dataset identity

The audited object is the GEARS/PertNet processed Norman file extracted from `norman.zip`.

| Field | Value |
| --- | --- |
| Local path | `data/norman/perturb_processed.h5ad` |
| SHA-256 | `23ffb0fac6a847ff927cf7509d80d85052bfefbfb97610786a2dafaaefa0b6a0` |
| Shape | 91,205 cells × 5,045 genes |
| `obs` | `condition`, `cell_type`, `dose_val`, `control`, `condition_name` |
| `var` | index `gene_id` (ENSG…); column `gene_name` |
| Main matrix `X` | CSR float32 log1p(per-cell scaled counts). `X` is **not** a stored extra layer. |
| Stored extra layer | `counts` only (integer-valued float32). Schema validation does **not** load the counts payload. |
| Feature order | **Stable `gene_id` values govern response-feature order.** Stored annotations are immutable. |

This report does **not** claim that the stored `counts` layer was reconstructed from the GEO UMI matrix.

---

## 3. Cell line: K562, high confidence

Independent sources agree the experiment used K562:

1. Norman et al., Science 2019: “K562 cells stably expressing the SunTag CRISPRa system”.
2. GEO series/sample: organism *Homo sapiens*; source name and characteristics name K562.
3. All 91,205 local barcodes join the GSE133344 raw identity table from that K562 experiment.

The only A549 evidence is the stored dummy column written by [Norman19 preprocessing cell 36](https://github.com/yhr91/GEARS_misc/blob/f88211870dfa89c38a2eedbd69ca1abd28a25f3c/data/preprocessing/Norman19.ipynb). [GEARS issue 95](https://github.com/snap-stanford/GEARS/issues/95#issuecomment-3216386250) confirms the dataloader does not use that attribute.

**Falsifier:** a different GEO accession whose barcodes match this file but whose samples are authenticated A549. No such join exists in the current evidence.

**Alternative that was checked and rejected:** “the file is A549 because `obs.cell_type` says so.” That label is literally assigned in preprocessor source, then echoed into `condition_name`.

---

## 4. Mechanism: SunTag CRISPRa, with a GEO KRAB conflict

The paper and the GEO series design describe **gain-of-function / overexpression / CRISPRa**. The paper names the **SunTag CRISPRa** reagent in K562.

GSM3906020 growth protocol nevertheless says “K562 cells expressing **dCas9-KRAB**”. That sentence conflicts with the paper and with the series overall design. Treating it as a GEO metadata error is a **strong inference**, not an author-confirmed correction. It does not outweigh the paper and series evidence for CRISPRa. Mechanism is **not** taken from that one line.

**Falsifier:** primary-author correction that this Perturb-seq plate actually used dCas9-KRAB. None is present in the downloaded sources.

---

## 5. Barcode groups (technical lanes)

Direct source `gemgroup` matches every local barcode suffix (`-1` … `-8`). Group counts:

| Group | Expression sample | Guide sample | Cells | Controls | Single-target cells | Double-target cells |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | GSM3906020 | GSM3906028 | 12409 | 1008 | 6548 | 4853 |
| 2 | GSM3906021 | GSM3906029 | 11422 | 928 | 6068 | 4426 |
| 3 | GSM3906022 | GSM3906030 | 11583 | 863 | 6216 | 4504 |
| 4 | GSM3906023 | GSM3906031 | 10629 | 874 | 5596 | 4159 |
| 5 | GSM3906024 | GSM3906032 | 10606 | 821 | 5635 | 4150 |
| 6 | GSM3906025 | GSM3906033 | 11152 | 976 | 5801 | 4375 |
| 7 | GSM3906026 | GSM3906034 | 11667 | 934 | 6308 | 4425 |
| 8 | GSM3906027 | GSM3906035 | 11737 | 949 | 6235 | 4553 |

Every group contains all 284 original condition labels and all 105 target genes. Full 105 × 8 counts: [target_gene_by_barcode_group.csv](target_gene_by_barcode_group.csv), [target_gene_by_barcode_group_wide.csv](target_gene_by_barcode_group_wide.csv). Lane summaries: [barcode_group_summary.csv](barcode_group_summary.csv).

Control-only Kruskal tests on retained-gene counts and detected features are significant but small (ε² ≈ 0.004). Control mean-profile pairwise correlations are 0.9983–0.9998 ([control_profile_correlations.csv](control_profile_correlations.csv)). These are descriptive technical-lane differences, not a donor effect and not a license to treat the eight groups as biological replicates.

GEO’s word “replicate” here means **lanes of one pooled experiment**.

---

## 6. Source QC

| Check | Result |
| --- | --- |
| Raw identity rows | 114,870 |
| Local barcodes in raw table | 91,205 / 91,205 |
| Condition match before RHOXF2 mapping | 90,443 |
| Condition match after documented RHOXF2 → RHOXF2BB | 91,205 |
| `gemgroup` match | 91,205 |
| `good_coverage` | all true |
| `number_of_cells` | all 1 |
| `cellranger_called` | 88,843 true / 2,362 false |
| Filtered-table matches | 88,843 |
| Filtered missing | 2,362 |
| Filtered guide conflicts | 5 ([filtered_source_conflicts.csv](filtered_source_conflicts.csv)) |
| Exact duplicate `X` rows | 0 (profile uniqueness only; does **not** prove absence of pairing) |
| Pairing IDs or before/after design | **Absent** from `obs`, `uns`, and joined source identity tables |

Zero duplicate `X` rows does **not** establish that cells lack before/after pairing. The conclusion that paired-cell counterfactuals are unavailable rests on the **absence of pairing IDs and pairing design**, not on row uniqueness.

The five filtered conflicts also have `good_coverage=false` and `number_of_cells=0` in the filtered table, so the filtered file is the wrong identity source for this H5AD.

Guide-assignment fields in the join (`read_count`, `UMI_count`, `coverage`) are source QC, not model inputs unless the user later approves them.

---

## 7. Naming anomalies

Keep stored tokens. Do not silently “fix” symbols.

| Stored / source | What is known | What is not known |
| --- | --- | --- |
| `RHOXF2` (source) / `RHOXF2BB` (stored) / `RHOXF2B` (author alias) | 762 cells; mapping is explicit in analyze_provenance, the GEARS `map_dict`, and the author notebook comment. `data_sharing` lists RHOXF2 and RHOXF2B with the **same** protospacer `GCATGCGTTGTCCTGTAGAC`. | Whether the molecular target is specifically RHOXF2, RHOXF2B, both, or a shared promoter. This is **dataset-local reconciliation**, not HGNC equivalence. |
| `C19orf26` / `CBARP`; `C3orf72` / `FOXL2NB`; `KIAA1804` / `RP5-862P8.2` | GEARS notebook remaps those Ensembl names **to** the stored tokens. | Whether later gene-nomenclature updates should replace the stored tokens. They must not be replaced without user approval. |
| `SLC38A2` | Assigned in the H5AD (280 single + 297 combinatorial cells; `ENSG00000134294`). | **No protospacer** in either author guide table. Left unresolved. |
| Control constructs | Three NegCtrl IDs collapse to stored `ctrl` (2,920 + 2,374 + 2,059 = 7,353). `NegCtrl1_…` was dropped by the preprocessor as the authors suggested. | Whether the three remaining controls are interchangeable reagents. |
| Suffix variants | `HOXC13`, `TGFBR2+IGDCC3`, and `ZBTB10` have `_1` / `_2` construct IDs that collapse to one condition each. | Whether those suffixes are distinct guides or sequencing duplicates. |

[verified_target_guide_sequences.csv](verified_target_guide_sequences.csv) records source gene, source guide ID, sequence, mapping rule, and unresolved status. [target_guide_sequences.csv](target_guide_sequences.csv) remains the old incomplete table.

`bulk_fitness` from the data_sharing table is an **outcome-derived fitness score**. It is external leakage if used as a covariate or feature.

---

## 8. Condition imbalance

Counts are from stored `condition` rows, not unused categorical levels. Sources: [perturbation_evidence.json](perturbation_evidence.json) (`count_summary_all`, `count_summary_noncontrol`, `canonical_count_summary_noncontrol`), [perturbation_counts.csv](perturbation_counts.csv), [canonical_perturbation_counts.csv](canonical_perturbation_counts.csv), and [DATA_AUDIT.md](DATA_AUDIT.md) §8.

| Aggregation | n | Min | Median | Max | Max/min | Extreme labels |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| All original conditions | 284 | 49 | 272 | 7,353 | 150.1 | max `ctrl` |
| Noncontrol original labels | 283 | 49 | 272 | 1,030 | 21.0 | min `C3orf72+FOXL2`; max `CEBPE+RUNX1T1` |
| Canonical noncontrol target sets | 236 | 49 | 294 | 1,641 | 33.5 | max `KLF1` (merges `KLF1+ctrl` and `ctrl+KLF1`) |

Cell-class totals: controls **7,353** / singles **48,407** / doubles **35,445**. Controls are **7,353 / 91,205 = 8.06%** of cells (ratio 1 : 11.40 vs 83,852 noncontrol cells). Cell representation is unequal; macro-averaging gives each condition equal weight rather than letting abundant conditions dominate.

---

## 9. Preprocessing

GEARS_misc `data/preprocessing/Norman19.ipynb` at commit `f88211870dfa89c38a2eedbd69ca1abd28a25f3c`:

1. Load `Norman2019_raw_all_genes.h5ad` (output: 93,658 × 24,665).
2. Drop `NegCtrl1_NegCtrl0__NegCtrl1_NegCtrl0` → 91,205 × 24,665.
3. Merge remaining NegCtrl constructs to `ctrl`; parse `guide_identity` with the second `__` half only.
4. Store a counts layer; `sc.pp.normalize_total(adata_new)` with default target; `sc.pp.log1p`; HVG `n_top_genes=5000` on all cells.
5. Force-keep perturbed gene names; remap CBARP/FOXL2NB/RP5-862P8.2/RHOXF2B onto stored tokens; subset → 91,205 × 5,045.
6. Write dummy `cell_type='A549'`, `dose_val`, `control`, and `uns` rankings.

scanpy 1.8.1 default `target_sum` is the median of nonzero per-cell totals. That median on the **24,665-gene** matrix is **unknown** here. A 10,000-total hypothesis on the retained 5,045 genes fails (only 1.96% of implied originals are near-integer). Full original-library reconstruction was **not** performed.

**HVG caveat:** 5,000 HVGs plus forced targets were selected on the full processed cohort. That is not a prospective gene-discovery claim. The same stored HVG feature set **may** be used in a retrospective prepared-feature benchmark if that caveat is stated. Existing `uns` DE rankings are **prohibited for fitting**.

---

## 10. Task A/B/C/D feasibility (proposed, not approved)

Seed `20260911`. Manifest: [proposed_split_assignments.csv](proposed_split_assignments.csv). Sizes from [provenance_evidence.json](provenance_evidence.json).

| Protocol | What it tests | Feasible? | Main limit |
| --- | --- | --- | --- |
| **A** random cell (within condition) | Interpolation / leakage-prone sanity | Yes, as a non-generalization check | Train/val/test share all 236 noncontrol target sets |
| **B** held-out combination | Predict unseen doubles from seen singles + other doubles | Yes | Does not hold out genes; singles of test-pair members remain in train |
| **C** held-out target gene | Unseen genes, including their combinations | Yes, after dropping 1,046 val–test bridge cells | Train genes may still appear as partners in val/test pairs |
| **D** held-out barcode group | Unseen technical lane | Yes, as a lane check | Not a biological-replicate or donor holdout; all perturbations remain visible |

---

## 11. PROPOSED primary protocol: B, combination held-out

**Proposal only. Not approved.**

All 105 single-target sets stay in train. 131 doubles are hashed into 26 test / 26 validation / 79 train. Controls are split 10% / 10% / remainder inside each barcode group.

| Split | Cells | Controls | Noncontrol cells | Noncontrol target sets | Singles | Doubles |
| --- | --- | --- | --- | --- | --- | --- |
| Train | 75,588 | 5,891 | 69,697 | 184 | 105 | 79 |
| Validation | 7,639 | 731 | 6,908 | 26 | 0 | 26 |
| Test | 7,978 | 731 | 7,247 | 26 | 0 | 26 |

Exact pair assignments: `B_pair_assignments` in [provenance_evidence.json](provenance_evidence.json).

---

## 12. PROPOSED prediction target and metrics

**Proposal only. Not approved.**

Proposed target: the **condition-level mean log-expression profile across all 5,045 measured genes**, expressed as a delta versus **training-split controls**. This is **not** a target-gene-only readout.

For canonical condition *C*:

```text
delta(C) = mean(X[cells in C, :]) - mean(X[training-split control cells, :])
```

`X` is the stored log1p per-cell-scaled matrix. Training-split controls are the only reference. Feature order is stored **`gene_id`**.

Paired-cell counterfactuals are unavailable because **pairing IDs and a before/after design are absent**. Zero duplicate `X` rows is only a uniqueness fact; it does not prove that pairing is absent. Do not invent paired twins.

Proposed (not approved) metrics if this target is later chosen:

- **Primary:** macro-condition **MSE** and **MAE** on the 5,045-gene delta profile.
- **Secondary:** Pearson correlation of predicted vs observed `delta` (**delta Pearson**).
- **Baseline skills:** **common-shift** and **additive**.

Do **not** use stored `uns` DE rankings or `bulk_fitness` as fitting inputs. Do not recompute test-expression-derived gene lists for training.

---

## 13. Leakage risks

| Risk | Why it matters |
| --- | --- |
| Full-data HVG selection | Not a prospective gene set; still usable as a retrospective prepared-feature set with that caveat |
| Existing `uns` DE lists | Prohibited for fitting; they encode full-cohort test-expression rankings |
| Random-cell split A as a “generalization” score | Same combinations appear in train and test |
| Using `bulk_fitness` | Outcome-derived external fitness leakage |
| Using filtered GEO identities | 2,362 missing cells and 5 wrong guides |
| Treating A549 as biology | Dummy preprocessor label |
| Treating eight groups as donors | They are 10x lanes of one pool |
| Counting `__` construct halves as four genes | Inflates perturbation cardinality |
| Collapsing RHOXF2/RHOXF2B/RHOXF2BB as proven molecular identity | Only a local label reconciliation |
| Synthesizing paired control cells | No pairing IDs or design; inventing twins would be fictional |
| Feeding test-set controls into the delta baseline | The proposed baseline is **training controls only** |

---

## 14. Unresolved scientific choices (user-reserved)

The schema and this report do **not** decide:

- Final task, split, metrics, model, or visualizations. The proposed target above is the full-profile condition mean delta, not a reserved alternative.
- Whether orientation labels (`GENE+ctrl` vs `ctrl+GENE`) are experimentally interchangeable.
- Whether to drop the 2,362 non-`cellranger_called` cells.
- Whether the three NegCtrl constructs should stay pooled.
- Biological interpretation of RHOXF2 / RHOXF2B / RHOXF2BB.
- Recovery of an `SLC38A2` guide sequence from another source.
- Whether a later protocol uses a different readout than the proposed full-profile delta.
- Whether to rebuild counts from GEO raw MTX (not done).

---

## 15. Falsification / alternative-explanation checks

| Hypothesis | Check already performed | What would falsify the current conclusion |
| --- | --- | --- |
| This H5AD is a different download | SHA-256 and shape match prior audits | A different digest after an untouched file |
| Barcodes are not Norman/GSE133344 | 91,205/91,205 raw join | Any local barcode absent from the raw table |
| Eight groups are donors | GEO text + `gemgroup` identity | Sample metadata naming eight donors, or gemgroup ≠ suffix |
| Cells are A549 | Dummy assignment + K562 paper/GEO + barcode join | Authenticated A549 source with these barcodes |
| Mechanism is CRISPRi because of “dCas9-KRAB” | Paper and series design contradict that sentence; KRAB-as-error is a strong inference | Author correction that KRAB was actually used |
| Filtered identities are safer | Filtered misses 2,362 and conflicts on 5 | A filtered table that matches all local barcodes and guides |
| RHOXF2 and RHOXF2B are different guides | Same protospacer in data_sharing | Distinct sequences plus a source that assigns both independently |
| Paired-cell deltas are available | Pairing IDs and before/after design are absent; 0 duplicate `X` rows is not that proof | A pairing ID, matched-cell design, or author-documented before/after pairing |

---

## 16. Reproducible validation (no training)

```powershell
python scripts/write_verified_guide_sequences.py
python scripts/validate_dataset_schema.py
python -m unittest discover -s tests -v
```

Do not rerun `scripts/download_norman.py` or execute downloaded upstream notebooks/code as project dependencies.
