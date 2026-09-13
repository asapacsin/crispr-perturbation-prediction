# AnnData data audit — Norman archive

Audited: 2026-09-11. **Completed; no model was trained.** The downloaded H5AD was read without changing it. Counts below describe this processed file, not necessarily the complete original experiment.

## Principal findings and evidence

| Requested finding | Conclusion | Direct evidence |
| --- | --- | --- |
| Physical observation | One cell-barcode-associated single-cell expression profile; interpreted as one captured cell. Singlet status is not independently established. | `obs.index.name == 'cell_barcode'`; 91,205 unique barcodes; one expression row and one condition assignment per barcode. Example: `AAACCTGAGGCATGTG-1`. |
| Variable | One measured gene-expression feature, identified by a gene ID and gene name. | `var.index.name == 'gene_id'`; 5,045 unique IDs and 5,045 unique `gene_name` values. Example: `ENSG00000187608` → `ISG15`. |
| Shape of X | **91,205 × 5,045**, cells × genes. | AnnData shape and HDF5 `X` shape both agree. CSR sparse matrix, `float32`. |
| Expression scale | **Natural-log1p-transformed, per-cell-scaled counts**, not raw counts or merely unlogged normalized counts. | All 37,317,477 nonzero entries were compared against the stored count layer; maximum absolute error in `X = log1p(count × scale[cell])` is 6.13492368e-07. |
| Important annotations | Five `obs` columns and one `var` column, detailed below. | Exact stored columns: `condition`, `cell_type`, `dose_val`, `control`, `condition_name`; `var['gene_name']`. None have missing values. |
| Perturbation representation | Gene-symbol tokens joined by `+`; `ctrl` marks a control component. | Observed labels include `TSC22D1+ctrl`, `ctrl+CEBPE`, `KLF1+MAP2K6`, and `ctrl`. |
| Controls | **7,353 cells** with `condition == 'ctrl'` and `control == 1`. | These predicates agree for every row. Their `dose_val` is `1` and `condition_name` is `A549_ctrl_1`. |
| Cells per perturbation | **284 original condition labels**, including control; **237 target sets**, including control, after explicit aggregation. | Full label counts and corresponding target-set totals appear in Appendix A; both original and aggregated counts sum to 91,205. |
| Batches | **Eight barcode suffix groups exist; experimental batch identities are not annotated.** | Barcode suffixes `-1` through `-8`; no batch/sample/library/donor/replicate column or such `uns` key. Group counts and limitations appear below. |
| Actual annotated targets | **105 distinct assigned target genes**, all measured in `var`. | Split observed condition labels on `+`, remove exactly `ctrl`, take distinct symbols, then exact-match against `var['gene_name']`: 105 matches, zero unmatched. Full list below. |

**Material identity conflict:** the supplied description calls this Norman et al. 2019, human K562, CRISPRa. However, **all 91,205 stored `cell_type` values are `A549`**, and all condition-name prefixes are `A549_`. This is a conflict between the supplied provenance and the file's annotations, not proof that the cells are biologically A549. The file contains no study citation, species, cell-line authentication, or CRISPR modality field that resolves it. Neither K562 identity nor activation versus knockout can be established from these annotations alone. No label was silently corrected.

## Source, integrity, and audit scope

- Requested local path: `data/norman/perturb_processed.h5ad`.
- Requested source: [Dataverse download, datafile 6154020](https://dataverse.harvard.edu/api/access/datafile/6154020).
- The response is a ZIP archive named `norman.zip`, not a bare H5AD. It was saved as `data/norman/norman.zip`; its size is **168,758,985 bytes**.
- Archive members, read from its directory: `norman/`, `norman/go.csv`, and `norman/perturb_processed.h5ad`. The H5AD member was extracted to the requested path. The archive's H5AD-member CRC32 is **a6b48d95**; reading the entire member through Python `zipfile` completed without a CRC error. Its uncompressed size is **2,228,610,012 bytes**. `go.csv` was not used to infer perturbation targets.
- Extracted H5AD SHA-256: `23ffb0fac6a847ff927cf7509d80d85052bfefbfb97610786a2dafaaefa0b6a0`. This records local file identity; no separately published SHA-256 was available for comparison.
- HDF5 top-level keys are exactly `X`, `layers`, `obs`, `uns`, `var`; the only stored layer is `counts`. There is no `raw`, `obsm`, `varm`, `obsp`, or processing-history object.
- Direct numerical scans covered the complete matrix and count layer. No fitting, prediction, dimensional reduction, differential-expression recomputation, or model training was performed.

## 1–2. Observations and variables

Each observation is a cell-barcode record with a transcript-expression vector and assigned perturbation. The index contains **91,205 unique, nonmissing barcodes**, each 18 characters long, including a numeric suffix. This is evidence for a single-cell profile as the unit of observation; it is not a perturbation average, a gene, or a biological replicate. The file provides no doublet score or singlet-validation field, so the physical claim is one *putative captured cell per barcode*, not a guarantee that every barcode arose from exactly one physical cell. Cell-line identity remains unresolved as noted above.

Each variable is an annotated gene feature. IDs have the `ENSG` naming pattern (for example `ENSG00000187608`), paired with gene names (for example `ISG15`). These are gene-expression features, not guide sequences: `var` contains only `gene_name`, and the matrix associates their expression with cell barcodes. Gene names include protein-coding-style symbols and noncoding-style names such as `LINC01342` and `RP11-34P13.8`; no biotype annotation is supplied. Both feature identifiers and names are unique. The file does not record the selection rule that produced these 5,045 features; it cannot be assumed that they are all genes originally measured or all highly variable genes.

Evidence: [complete observation annotations](obs_metadata.csv), [complete feature annotations](var_metadata.csv), [structural and numerical evidence](audit_evidence.json).

## 3–4. Matrix shape, counts, and transformation

`X` has **460,129,225 possible entries**, of which **37,317,477 are nonzero** (8.1102% density; 91.8898% zeros). It has no empty rows or empty columns. All values are finite and nonnegative. Implicit sparse zeros are included in these conclusions.

| Quantity | `X` | `layers['counts']` |
| --- | --- | --- |
| Shape | 91,205 × 5,045 | 91,205 × 5,045 |
| Storage read into memory | CSR sparse, float32 | Dense NumPy array, float32 |
| Nonzero entries | 37,317,477 | 37,317,477 |
| Minimum including zeros | 0 | 0 |
| Maximum | 8.9045801163 | 3,718 |
| Entries more than 10⁻⁶ from an integer | 37,317,464 | 0 |
| Negative or nonfinite entries | 0 | 0 |
| Per-row sum: min / median / max | 289.274994 / 538.095398 / 930.651550 | 322 / 2,999 / 11,831 |

The layer named `counts` contains **untransformed, nonnegative integer-valued counts**, despite its float32 dtype. It is the stored raw-count-like representation. Interpreting these specifically as raw UMI molecule counts is consistent with the supplied Perturb-seq provenance, but the file itself contains no UMI/read-count provenance or correction history. Integer values and a layer name alone do not prove the original molecule-counting procedure.

The transformation classification is supported numerically, not just by fractional values. For each cell, the audit derived a scale from its first nonzero entry as `expm1(X[i,j]) / counts[i,j]`, then checked every other nonzero entry in that cell against:

```text
X[i,j] = ln(1 + scale[i] × counts[i,j])
```

All **37,317,477** nonzero entries match within **2 × 10⁻⁶** absolute error; the worst error is **6.13492368e-07**. Zero locations match exactly. The per-cell scales range from **0.2474850895 to 6.5645891100**, with median **1.0000000355**.

For the first barcode, `AAACCTGAGGCATGTG-1`, directly observed entries are:

| Gene | Stored count | X | `expm1(X)/count` |
| --- | --- | --- | --- |
| ISG15 | 1 | 0.6815017462 | 0.9768442225 |
| ENO1 | 4 | 1.5907396078 | 0.9768442802 |
| SRM | 2 | 1.0830547810 | 0.9768443278 |

Thus the data support **per-cell scaling followed by natural log1p**. They do not establish a particular normalization function, denominator gene set, or target total. The retained genes' `expm1(X)` row sums range from **1,964.0674 to 9,479.6816** (median **2,963.7360**), so those sums are not a constant 10,000. No `uns['log1p']` or normalization settings are stored. A tested 10,000-total hypothesis does not recover consistently integer original totals (only **1.957%** of `10000/scale` values lie within 0.01 of an integer); this does not reconstruct or prove the original normalization recipe.

There is no AnnData `raw` snapshot. The installed AnnData version exposes `layers[None]` as an API alias for `X`; HDF5 inspection confirms that this is **not** a second stored layer named `null` or `None`.

Evidence: [audit_evidence.json](audit_evidence.json), [full-matrix transformation validation](validation_evidence.json).

## 5. Complete annotation inventory

| Location | Field | Stored/read dtype | Observed values and interpretation | Missing |
| --- | --- | --- | --- | --- |
| `obs.index` | `cell_barcode` | string identifiers | 91,205 unique cell barcodes; preserve numeric suffixes. | 0 |
| `obs` | `condition` | categorical | 284 observed labels, the primary perturbation assignment. Examples: `ctrl`, `TSC22D1+ctrl`, `ctrl+CEBPE`, `KLF1+MAP2K6`. | 0 |
| `obs` | `cell_type` | categorical | One value: `A549` in 91,205 cells; conflicts with supplied K562 description. | 0 |
| `obs` | `dose_val` | categorical | `1` in 7,353 control cells; `1+1` in all 83,852 other cells. No dose units supplied; not evidence of measured physical dose or perturbation efficacy. | 0 |
| `obs` | `control` | int64 | Binary: 1 in 7,353 cells; 0 in 83,852 cells. Exactly agrees with `condition == 'ctrl'`. | 0 |
| `obs` | `condition_name` | categorical | 284 labels. Exactly equals `cell_type + '_' + condition + '_' + dose_val` for every row. | 0 |
| `var.index` | `gene_id` | string identifiers | 5,045 unique `ENSG...` identifiers. | 0 |
| `var` | `gene_name` | categorical | 5,045 unique names; exact mapping to the feature index. | 0 |

These are **all** columns, not a selected subset. There are no guide IDs, guide sequences, guide-UMI counts, assignment-confidence scores, sample IDs, batch IDs, library IDs, donor IDs, replicate IDs, mitochondrial fractions, total-count annotations, gene biotypes, or highly-variable flags among them. The expression matrix can support derived summaries, but no such missing provenance fields were invented.

The five `uns` keys hold previously computed per-condition arrays:

| `uns` key | Observed structure | Interpretation supported by the stored content |
| --- | --- | --- |
| `non_dropout_gene_idx` | 283 condition keys; integer arrays of length 2,731–4,063; values 0–5,044. | Stored gene-position subsets, labeled non-dropout. Exact selection algorithm is not recorded. |
| `non_zeros_gene_idx` | 283 condition keys; integer arrays of length 1,411–2,863; values 0–5,044. | Stored gene-position subsets, labeled nonzero. Exact selection algorithm is not recorded. |
| `rank_genes_groups_cov_all` | 283 condition keys; 5,045 gene IDs per key. All IDs occur in `var.index`. | Stored gene-ID orderings, named as rankings; no associated test settings, scores, or provenance are present. |
| `top_non_dropout_de_20` | 283 condition keys; 20 gene IDs per key, all present in `var.index`. | Stored top-20 lists named as non-dropout differential-expression results. |
| `top_non_zero_de_20` | 283 condition keys; 20 gene IDs per key, all present in `var.index`. | Stored top-20 lists named as nonzero differential-expression results. |

Every `uns` dictionary's keys match the noncontrol `condition_name` values exactly; the only condition absent is `A549_ctrl_1`. Example key: `A549_AHR+FEV_1+1`. These arrays are not perturbation-target lists. Their biological/statistical validity was not assumed or re-estimated. Evidence and representative array values: [validation_evidence.json](validation_evidence.json).

## 6–7. Perturbations and controls

A stored `condition` is parsed literally using `+`. The token `ctrl` is removed only when deriving the set of target genes:

| Example original label | Interpretation from annotation | Distinct gene targets |
| --- | --- | --- |
| `ctrl` | Control-labeled cell | none |
| `TSC22D1+ctrl` | Single-target condition with a control component | TSC22D1 |
| `ctrl+CEBPE` | Single-target condition with the control component first | CEBPE |
| `KLF1+MAP2K6` | Two-target condition | KLF1, MAP2K6 |

The annotation supports zero, one, or two distinct noncontrol gene tokens per cell. It does not specify guide sequences, whether the control guide is non-targeting versus another control design, or successful activation of any assigned target. In particular, a cell with `GENE+ctrl` or `ctrl+GENE` is **not** a control cell: every such row has `control == 0`.

There are **152 original single-target labels**, representing **105 target genes**, because **47 genes** have both `GENE+ctrl` and `ctrl+GENE` labels. The original orientation is preserved in the label-level counts. Aggregated target sets ignore order and remove `ctrl`; this is an explicit descriptive aggregation, not evidence that different construct orientations are experimentally interchangeable. Example: `KLF1+ctrl` has **997** cells and `ctrl+KLF1` has **644**, totaling **1,641 KLF1-targeted cells**. All **131 two-target labels** map to different unordered two-gene sets; no two-target labels collapse together in this file.

Control checks across all observations show exact agreement of `condition == 'ctrl'`, `control == 1`, and `dose_val == '1'`. All other rows have `dose_val == '1+1'`. Evidence: [perturbation_evidence.json](perturbation_evidence.json), [original metadata](obs_metadata.csv).

## 8. Cells per perturbation

| Class | Original labels | Distinct target sets | Cells |
| --- | --- | --- | --- |
| Control | 1 | 1 | 7,353 |
| Single target | 152 | 105 | 48,407 |
| Two targets | 131 | 131 | 35,445 |
| Total | 284 | 237 | 91,205 |

Noncontrol original-label sizes range from **49 to 1,030 cells**, with median **272** and mean **296.297**. The smallest is `C3orf72+FOXL2` (49), followed by `CBL+UBASH3A` (50); the largest is `CEBPE+RUNX1T1` (1,030). After target-set aggregation, the largest noncontrol set is `KLF1` (1,641); the median noncontrol set size is **294**. These counts establish that representation is unequal across conditions.

**Every original label and count is included in Appendix A**, together with its aggregated target-set count. Machine-readable tables: [284 original-label counts](perturbation_counts.csv), [237 target-set counts with original-label mapping](canonical_perturbation_counts.csv). Counts are based on observed rows, not unused categorical levels. A pair-target cell counts once in these tables; its two targets are counted separately only in the per-gene exposure table.

## 9. Batch evidence and limits

There is no explicit experimental batch variable. However, splitting the barcode at its final hyphen yields exactly eight groups:

| Barcode suffix | Cells | Controls | Observed condition labels |
| --- | --- | --- | --- |
| 1 | 12409 | 1008 | 284 |
| 2 | 11422 | 928 | 284 |
| 3 | 11583 | 863 | 284 |
| 4 | 10629 | 874 | 284 |
| 5 | 10606 | 821 | 284 |
| 6 | 11152 | 976 | 284 |
| 7 | 11667 | 934 | 284 |
| 8 | 11737 | 949 | 284 |

All **284 condition labels occur in every suffix group**, and controls occur in every group. This is evidence of eight barcode partitions and is consistent with multiple source libraries or groups, but the file does not map suffixes to collection dates, lanes, libraries, or biological/technical replicates. Therefore **“no batches exist” is unsupported**, and **“eight experimentally verified batches” would also overstate the evidence**. The suffix is a candidate technical grouping, not an authenticated batch annotation. The audit did not test expression batch effects.

The complete barcodes are unique. Removing suffixes introduces **4,768 duplicate occurrences beyond the first**, so suffixes must be preserved as part of cell identity. Evidence: [suffix summary](barcode_suffix_counts.csv), [condition × suffix counts](condition_by_barcode_suffix.csv), [validation_evidence.json](validation_evidence.json).

## 10. Genes actually targeted according to the assignments

The **105 distinct symbols below** come exclusively from observed perturbation labels, after removing `ctrl`. All have at least one single-target condition in this file and all match `var['gene_name']` exactly. The other **4,940 measured features** are not named as targets in any observed condition. A measured gene or a gene in a stored ranking was not counted as a target on that basis.

`AHR`, `ARID1A`, `ARRDC3`, `ATL1`, `BAK1`, `BCL2L11`, `BCORL1`, `BPGM`, `C19orf26`, `C3orf72`, `CBFA2T3`, `CBL`, `CDKN1A`, `CDKN1B`, `CDKN1C`.

`CEBPA`, `CEBPB`, `CEBPE`, `CELF2`, `CITED1`, `CKS1B`, `CLDN6`, `CNN1`, `CNNM4`, `COL1A1`, `COL2A1`, `CSRNP1`, `DLX2`, `DUSP9`, `EGR1`.

`ELMSAN1`, `ETS2`, `FEV`, `FOSB`, `FOXA1`, `FOXA3`, `FOXF1`, `FOXL2`, `FOXO4`, `GLB1L2`, `HES7`, `HK2`, `HNF4A`, `HOXA13`, `HOXB9`.

`HOXC13`, `IER5L`, `IGDCC3`, `IKZF3`, `IRF1`, `ISL2`, `JUN`, `KIAA1804`, `KIF18B`, `KIF2C`, `KLF1`, `KMT2A`, `LHX1`, `LYL1`, `MAML2`.

`MAP2K3`, `MAP2K6`, `MAP4K3`, `MAP4K5`, `MAP7D1`, `MAPK1`, `MEIS1`, `MIDN`, `NCL`, `NIT1`, `OSR2`, `PLK4`, `POU3F2`, `PRDM1`, `PRTG`.

`PTPN1`, `PTPN12`, `PTPN13`, `PTPN9`, `RHOXF2BB`, `RREB1`, `RUNX1T1`, `S1PR2`, `SAMD1`, `SET`, `SGK1`, `SLC38A2`, `SLC4A1`, `SLC6A9`, `SNAI1`.

`SPI1`, `STIL`, `TBX2`, `TBX3`, `TGFBR2`, `TMSB4X`, `TP73`, `TSC22D1`, `UBASH3A`, `UBASH3B`, `ZBTB1`, `ZBTB10`, `ZBTB25`, `ZC3HAV1`, `ZNF318`.

The symbols are retained as stored, including names such as `C19orf26`, `KIAA1804`, and `RHOXF2BB`; no external alias conversion was guessed. [target_genes.csv](target_genes.csv) gives each target's exact feature ID, number of cells in single-target conditions, number in two-target conditions, total assigned cells, and number of original labels. Per-gene totals overlap for two-target cells and must not be summed as a count of unique cells.

“Targeted” here means **assigned an annotated perturbation against that gene**. Successful molecular activation, guide identity, and off-target effects cannot be established from the supplied fields.

## Reproducibility and evidence files

The audit can be rerun from the project root, with the H5AD already present:

```powershell
python scripts/inspect_norman.py
python scripts/audit_perturbations.py
python scripts/validate_norman.py
python scripts/write_data_audit.py
```

Inspected with `anndata 0.13.3.post0`, `numpy 2.5.2`, `scipy 1.18.1`, `pandas 3.0.5`, `h5py 3.16.0`.

The scripts inspect stored data and perform deterministic summaries and numerical identity checks; they do not train models. Numerical evidence is retained in [audit_evidence.json](audit_evidence.json), [perturbation_evidence.json](perturbation_evidence.json), and [validation_evidence.json](validation_evidence.json). Complete `obs` and `var` exports and all count tables are linked above. The download archive's member manifest is retained in `data/norman/archive_manifest.json`.

## Appendix A. Every original perturbation label and cell count

The target-set total is repeated when two original labels map to the same target set. **Sum the original-label cells column to obtain 91,205; do not sum the repeated target-set totals.** Use the separate canonical table for unique target-set totals.

| Original condition | Cells | Target set | Target-set cells |
| --- | --- | --- | --- |
| AHR+FEV | 264 | AHR+FEV | 264 |
| AHR+KLF1 | 412 | AHR+KLF1 | 412 |
| AHR+ctrl | 479 | AHR | 479 |
| ARID1A+ctrl | 182 | ARID1A | 182 |
| ARRDC3+ctrl | 405 | ARRDC3 | 405 |
| ATL1+ctrl | 305 | ATL1 | 305 |
| BAK1+ctrl | 534 | BAK1 | 1171 |
| BCL2L11+BAK1 | 153 | BAK1+BCL2L11 | 153 |
| BCL2L11+TGFBR2 | 382 | BCL2L11+TGFBR2 | 382 |
| BCL2L11+ctrl | 463 | BCL2L11 | 463 |
| BCORL1+ctrl | 456 | BCORL1 | 456 |
| BPGM+SAMD1 | 240 | BPGM+SAMD1 | 240 |
| BPGM+ZBTB1 | 283 | BPGM+ZBTB1 | 283 |
| BPGM+ctrl | 393 | BPGM | 393 |
| C19orf26+ctrl | 480 | C19orf26 | 750 |
| C3orf72+FOXL2 | 49 | C3orf72+FOXL2 | 49 |
| C3orf72+ctrl | 217 | C3orf72 | 217 |
| CBFA2T3+ctrl | 288 | CBFA2T3 | 531 |
| CBL+CNN1 | 288 | CBL+CNN1 | 288 |
| CBL+PTPN12 | 257 | CBL+PTPN12 | 257 |
| CBL+PTPN9 | 234 | CBL+PTPN9 | 234 |
| CBL+TGFBR2 | 156 | CBL+TGFBR2 | 156 |
| CBL+UBASH3A | 50 | CBL+UBASH3A | 50 |
| CBL+UBASH3B | 326 | CBL+UBASH3B | 326 |
| CBL+ctrl | 538 | CBL | 538 |
| CDKN1A+ctrl | 144 | CDKN1A | 275 |
| CDKN1B+CDKN1A | 98 | CDKN1A+CDKN1B | 98 |
| CDKN1B+ctrl | 268 | CDKN1B | 413 |
| CDKN1C+CDKN1A | 80 | CDKN1A+CDKN1C | 80 |
| CDKN1C+CDKN1B | 79 | CDKN1B+CDKN1C | 79 |
| CDKN1C+ctrl | 165 | CDKN1C | 165 |
| CEBPA+ctrl | 322 | CEBPA | 580 |
| CEBPB+CEBPA | 75 | CEBPA+CEBPB | 75 |
| CEBPB+MAPK1 | 337 | CEBPB+MAPK1 | 337 |
| CEBPB+OSR2 | 188 | CEBPB+OSR2 | 188 |
| CEBPB+PTPN12 | 266 | CEBPB+PTPN12 | 266 |
| CEBPB+ctrl | 240 | CEBPB | 422 |
| CEBPE+CEBPA | 240 | CEBPA+CEBPE | 240 |
| CEBPE+CEBPB | 111 | CEBPB+CEBPE | 111 |
| CEBPE+CNN1 | 194 | CEBPE+CNN1 | 194 |
| CEBPE+KLF1 | 389 | CEBPE+KLF1 | 389 |
| CEBPE+PTPN12 | 317 | CEBPE+PTPN12 | 317 |
| CEBPE+RUNX1T1 | 1030 | CEBPE+RUNX1T1 | 1030 |
| CEBPE+SPI1 | 175 | CEBPE+SPI1 | 175 |
| CEBPE+ctrl | 473 | CEBPE | 1039 |
| CELF2+ctrl | 388 | CELF2 | 388 |
| CITED1+ctrl | 169 | CITED1 | 169 |
| CKS1B+ctrl | 189 | CKS1B | 189 |
| CLDN6+ctrl | 281 | CLDN6 | 544 |
| CNN1+MAPK1 | 325 | CNN1+MAPK1 | 325 |
| CNN1+UBASH3A | 144 | CNN1+UBASH3A | 144 |
| CNN1+ctrl | 236 | CNN1 | 636 |
| CNNM4+ctrl | 376 | CNNM4 | 376 |
| COL1A1+ctrl | 257 | COL1A1 | 257 |
| COL2A1+ctrl | 247 | COL2A1 | 474 |
| CSRNP1+ctrl | 354 | CSRNP1 | 354 |
| DLX2+ctrl | 316 | DLX2 | 606 |
| DUSP9+ETS2 | 698 | DUSP9+ETS2 | 698 |
| DUSP9+IGDCC3 | 323 | DUSP9+IGDCC3 | 323 |
| DUSP9+KLF1 | 320 | DUSP9+KLF1 | 320 |
| DUSP9+MAPK1 | 247 | DUSP9+MAPK1 | 247 |
| DUSP9+PRTG | 238 | DUSP9+PRTG | 238 |
| DUSP9+SNAI1 | 168 | DUSP9+SNAI1 | 168 |
| DUSP9+ctrl | 662 | DUSP9 | 662 |
| EGR1+ctrl | 261 | EGR1 | 261 |
| ELMSAN1+ctrl | 353 | ELMSAN1 | 783 |
| ETS2+CEBPE | 320 | CEBPE+ETS2 | 320 |
| ETS2+CNN1 | 785 | CNN1+ETS2 | 785 |
| ETS2+IGDCC3 | 365 | ETS2+IGDCC3 | 365 |
| ETS2+IKZF3 | 388 | ETS2+IKZF3 | 388 |
| ETS2+MAP7D1 | 265 | ETS2+MAP7D1 | 265 |
| ETS2+MAPK1 | 432 | ETS2+MAPK1 | 432 |
| ETS2+PRTG | 354 | ETS2+PRTG | 354 |
| ETS2+ctrl | 375 | ETS2 | 1031 |
| FEV+CBFA2T3 | 153 | CBFA2T3+FEV | 153 |
| FEV+ISL2 | 272 | FEV+ISL2 | 272 |
| FEV+MAP7D1 | 246 | FEV+MAP7D1 | 246 |
| FEV+ctrl | 173 | FEV | 647 |
| FOSB+CEBPB | 85 | CEBPB+FOSB | 85 |
| FOSB+CEBPE | 142 | CEBPE+FOSB | 142 |
| FOSB+IKZF3 | 204 | FOSB+IKZF3 | 204 |
| FOSB+OSR2 | 324 | FOSB+OSR2 | 324 |
| FOSB+PTPN12 | 291 | FOSB+PTPN12 | 291 |
| FOSB+UBASH3B | 363 | FOSB+UBASH3B | 363 |
| FOSB+ctrl | 502 | FOSB | 502 |
| FOXA1+FOXF1 | 237 | FOXA1+FOXF1 | 237 |
| FOXA1+FOXL2 | 179 | FOXA1+FOXL2 | 179 |
| FOXA1+HOXB9 | 266 | FOXA1+HOXB9 | 266 |
| FOXA1+ctrl | 426 | FOXA1 | 761 |
| FOXA3+FOXA1 | 181 | FOXA1+FOXA3 | 181 |
| FOXA3+FOXF1 | 146 | FOXA3+FOXF1 | 146 |
| FOXA3+FOXL2 | 113 | FOXA3+FOXL2 | 113 |
| FOXA3+HOXB9 | 301 | FOXA3+HOXB9 | 301 |
| FOXA3+ctrl | 409 | FOXA3 | 409 |
| FOXF1+FOXL2 | 173 | FOXF1+FOXL2 | 173 |
| FOXF1+HOXB9 | 173 | FOXF1+HOXB9 | 173 |
| FOXF1+ctrl | 448 | FOXF1 | 715 |
| FOXL2+HOXB9 | 87 | FOXL2+HOXB9 | 87 |
| FOXL2+MEIS1 | 177 | FOXL2+MEIS1 | 177 |
| FOXL2+ctrl | 303 | FOXL2 | 544 |
| FOXO4+ctrl | 189 | FOXO4 | 189 |
| GLB1L2+ctrl | 378 | GLB1L2 | 378 |
| HES7+ctrl | 122 | HES7 | 122 |
| HK2+ctrl | 249 | HK2 | 249 |
| HNF4A+ctrl | 202 | HNF4A | 202 |
| HOXA13+ctrl | 230 | HOXA13 | 230 |
| HOXB9+ctrl | 356 | HOXB9 | 612 |
| HOXC13+ctrl | 151 | HOXC13 | 343 |
| IER5L+ctrl | 132 | IER5L | 327 |
| IGDCC3+MAPK1 | 293 | IGDCC3+MAPK1 | 293 |
| IGDCC3+PRTG | 108 | IGDCC3+PRTG | 108 |
| IGDCC3+ZBTB25 | 99 | IGDCC3+ZBTB25 | 99 |
| IGDCC3+ctrl | 284 | IGDCC3 | 507 |
| IKZF3+ctrl | 374 | IKZF3 | 583 |
| IRF1+SET | 367 | IRF1+SET | 367 |
| IRF1+ctrl | 331 | IRF1 | 331 |
| ISL2+ctrl | 272 | ISL2 | 485 |
| JUN+CEBPA | 70 | CEBPA+JUN | 70 |
| JUN+CEBPB | 52 | CEBPB+JUN | 52 |
| JUN+ctrl | 235 | JUN | 235 |
| KIAA1804+ctrl | 212 | KIAA1804 | 212 |
| KIF18B+KIF2C | 83 | KIF18B+KIF2C | 83 |
| KIF18B+ctrl | 212 | KIF18B | 212 |
| KIF2C+ctrl | 270 | KIF2C | 518 |
| KLF1+BAK1 | 323 | BAK1+KLF1 | 323 |
| KLF1+CEBPA | 308 | CEBPA+KLF1 | 308 |
| KLF1+CLDN6 | 215 | CLDN6+KLF1 | 215 |
| KLF1+COL2A1 | 251 | COL2A1+KLF1 | 251 |
| KLF1+FOXA1 | 323 | FOXA1+KLF1 | 323 |
| KLF1+MAP2K6 | 222 | KLF1+MAP2K6 | 222 |
| KLF1+TGFBR2 | 263 | KLF1+TGFBR2 | 263 |
| KLF1+ctrl | 997 | KLF1 | 1641 |
| KMT2A+ctrl | 185 | KMT2A | 185 |
| LHX1+ELMSAN1 | 547 | ELMSAN1+LHX1 | 547 |
| LHX1+ctrl | 342 | LHX1 | 342 |
| LYL1+CEBPB | 159 | CEBPB+LYL1 | 159 |
| LYL1+IER5L | 547 | IER5L+LYL1 | 547 |
| LYL1+ctrl | 360 | LYL1 | 360 |
| MAML2+ctrl | 508 | MAML2 | 508 |
| MAP2K3+ELMSAN1 | 484 | ELMSAN1+MAP2K3 | 484 |
| MAP2K3+IKZF3 | 220 | IKZF3+MAP2K3 | 220 |
| MAP2K3+MAP2K6 | 458 | MAP2K3+MAP2K6 | 458 |
| MAP2K3+SLC38A2 | 297 | MAP2K3+SLC38A2 | 297 |
| MAP2K3+ctrl | 458 | MAP2K3 | 458 |
| MAP2K6+ELMSAN1 | 365 | ELMSAN1+MAP2K6 | 365 |
| MAP2K6+IKZF3 | 251 | IKZF3+MAP2K6 | 251 |
| MAP2K6+SPI1 | 242 | MAP2K6+SPI1 | 242 |
| MAP2K6+ctrl | 363 | MAP2K6 | 741 |
| MAP4K3+ctrl | 299 | MAP4K3 | 299 |
| MAP4K5+ctrl | 300 | MAP4K5 | 300 |
| MAP7D1+ctrl | 277 | MAP7D1 | 635 |
| MAPK1+IKZF3 | 264 | IKZF3+MAPK1 | 264 |
| MAPK1+PRTG | 380 | MAPK1+PRTG | 380 |
| MAPK1+TGFBR2 | 389 | MAPK1+TGFBR2 | 389 |
| MAPK1+ctrl | 279 | MAPK1 | 660 |
| MEIS1+ctrl | 292 | MEIS1 | 682 |
| MIDN+ctrl | 256 | MIDN | 256 |
| NCL+ctrl | 197 | NCL | 197 |
| NIT1+ctrl | 165 | NIT1 | 165 |
| OSR2+ctrl | 407 | OSR2 | 858 |
| PLK4+STIL | 66 | PLK4+STIL | 66 |
| PLK4+ctrl | 97 | PLK4 | 97 |
| POU3F2+CBFA2T3 | 202 | CBFA2T3+POU3F2 | 202 |
| POU3F2+FOXL2 | 286 | FOXL2+POU3F2 | 286 |
| POU3F2+ctrl | 312 | POU3F2 | 312 |
| PRDM1+CBFA2T3 | 109 | CBFA2T3+PRDM1 | 109 |
| PRDM1+ctrl | 273 | PRDM1 | 273 |
| PRTG+ctrl | 285 | PRTG | 581 |
| PTPN1+ctrl | 342 | PTPN1 | 342 |
| PTPN12+OSR2 | 283 | OSR2+PTPN12 | 283 |
| PTPN12+PTPN9 | 295 | PTPN12+PTPN9 | 295 |
| PTPN12+SNAI1 | 244 | PTPN12+SNAI1 | 244 |
| PTPN12+UBASH3A | 255 | PTPN12+UBASH3A | 255 |
| PTPN12+ZBTB25 | 257 | PTPN12+ZBTB25 | 257 |
| PTPN12+ctrl | 194 | PTPN12 | 358 |
| PTPN13+ctrl | 226 | PTPN13 | 226 |
| PTPN9+ctrl | 164 | PTPN9 | 363 |
| RHOXF2BB+SET | 267 | RHOXF2BB+SET | 267 |
| RHOXF2BB+ZBTB25 | 205 | RHOXF2BB+ZBTB25 | 205 |
| RHOXF2BB+ctrl | 290 | RHOXF2BB | 290 |
| RREB1+ctrl | 281 | RREB1 | 281 |
| RUNX1T1+ctrl | 275 | RUNX1T1 | 631 |
| S1PR2+ctrl | 446 | S1PR2 | 446 |
| SAMD1+PTPN12 | 292 | PTPN12+SAMD1 | 292 |
| SAMD1+TGFBR2 | 59 | SAMD1+TGFBR2 | 59 |
| SAMD1+UBASH3B | 209 | SAMD1+UBASH3B | 209 |
| SAMD1+ZBTB1 | 164 | SAMD1+ZBTB1 | 164 |
| SAMD1+ctrl | 252 | SAMD1 | 430 |
| SET+CEBPE | 529 | CEBPE+SET | 529 |
| SET+KLF1 | 541 | KLF1+SET | 541 |
| SET+ctrl | 512 | SET | 790 |
| SGK1+S1PR2 | 312 | S1PR2+SGK1 | 312 |
| SGK1+TBX2 | 287 | SGK1+TBX2 | 287 |
| SGK1+TBX3 | 243 | SGK1+TBX3 | 243 |
| SGK1+ctrl | 446 | SGK1 | 446 |
| SLC4A1+ctrl | 853 | SLC4A1 | 853 |
| SLC6A9+ctrl | 441 | SLC6A9 | 441 |
| SNAI1+DLX2 | 107 | DLX2+SNAI1 | 107 |
| SNAI1+UBASH3B | 256 | SNAI1+UBASH3B | 256 |
| SNAI1+ctrl | 220 | SNAI1 | 380 |
| SPI1+ctrl | 170 | SPI1 | 506 |
| STIL+ctrl | 127 | STIL | 268 |
| TBX2+ctrl | 269 | TBX2 | 564 |
| TBX3+TBX2 | 969 | TBX2+TBX3 | 969 |
| TBX3+ctrl | 289 | TBX3 | 559 |
| TGFBR2+C19orf26 | 274 | C19orf26+TGFBR2 | 274 |
| TGFBR2+ETS2 | 275 | ETS2+TGFBR2 | 275 |
| TGFBR2+IGDCC3 | 262 | IGDCC3+TGFBR2 | 262 |
| TGFBR2+PRTG | 208 | PRTG+TGFBR2 | 208 |
| TGFBR2+ctrl | 209 | TGFBR2 | 537 |
| TMSB4X+BAK1 | 267 | BAK1+TMSB4X | 267 |
| TMSB4X+ctrl | 462 | TMSB4X | 462 |
| TP73+ctrl | 218 | TP73 | 218 |
| TSC22D1+ctrl | 376 | TSC22D1 | 376 |
| UBASH3A+ctrl | 371 | UBASH3A | 686 |
| UBASH3B+CNN1 | 355 | CNN1+UBASH3B | 355 |
| UBASH3B+OSR2 | 677 | OSR2+UBASH3B | 677 |
| UBASH3B+PTPN12 | 260 | PTPN12+UBASH3B | 260 |
| UBASH3B+PTPN9 | 282 | PTPN9+UBASH3B | 282 |
| UBASH3B+UBASH3A | 426 | UBASH3A+UBASH3B | 426 |
| UBASH3B+ZBTB25 | 224 | UBASH3B+ZBTB25 | 224 |
| UBASH3B+ctrl | 470 | UBASH3B | 983 |
| ZBTB1+ctrl | 315 | ZBTB1 | 602 |
| ZBTB10+DLX2 | 74 | DLX2+ZBTB10 | 74 |
| ZBTB10+ELMSAN1 | 161 | ELMSAN1+ZBTB10 | 161 |
| ZBTB10+PTPN12 | 230 | PTPN12+ZBTB10 | 230 |
| ZBTB10+SNAI1 | 85 | SNAI1+ZBTB10 | 85 |
| ZBTB10+ctrl | 145 | ZBTB10 | 145 |
| ZBTB25+ctrl | 343 | ZBTB25 | 590 |
| ZC3HAV1+CEBPA | 124 | CEBPA+ZC3HAV1 | 124 |
| ZC3HAV1+CEBPE | 346 | CEBPE+ZC3HAV1 | 346 |
| ZC3HAV1+HOXC13 | 543 | HOXC13+ZC3HAV1 | 543 |
| ZC3HAV1+ctrl | 436 | ZC3HAV1 | 436 |
| ZNF318+FOXL2 | 197 | FOXL2+ZNF318 | 197 |
| ZNF318+ctrl | 541 | ZNF318 | 541 |
| ctrl | 7353 | ctrl | 7353 |
| ctrl+BAK1 | 637 | BAK1 | 1171 |
| ctrl+C19orf26 | 270 | C19orf26 | 750 |
| ctrl+CBFA2T3 | 243 | CBFA2T3 | 531 |
| ctrl+CDKN1A | 131 | CDKN1A | 275 |
| ctrl+CDKN1B | 145 | CDKN1B | 413 |
| ctrl+CEBPA | 258 | CEBPA | 580 |
| ctrl+CEBPB | 182 | CEBPB | 422 |
| ctrl+CEBPE | 566 | CEBPE | 1039 |
| ctrl+CLDN6 | 263 | CLDN6 | 544 |
| ctrl+CNN1 | 400 | CNN1 | 636 |
| ctrl+COL2A1 | 227 | COL2A1 | 474 |
| ctrl+DLX2 | 290 | DLX2 | 606 |
| ctrl+ELMSAN1 | 430 | ELMSAN1 | 783 |
| ctrl+ETS2 | 656 | ETS2 | 1031 |
| ctrl+FEV | 474 | FEV | 647 |
| ctrl+FOXA1 | 335 | FOXA1 | 761 |
| ctrl+FOXF1 | 267 | FOXF1 | 715 |
| ctrl+FOXL2 | 241 | FOXL2 | 544 |
| ctrl+HOXB9 | 256 | HOXB9 | 612 |
| ctrl+HOXC13 | 192 | HOXC13 | 343 |
| ctrl+IER5L | 195 | IER5L | 327 |
| ctrl+IGDCC3 | 223 | IGDCC3 | 507 |
| ctrl+IKZF3 | 209 | IKZF3 | 583 |
| ctrl+ISL2 | 213 | ISL2 | 485 |
| ctrl+KIF2C | 248 | KIF2C | 518 |
| ctrl+KLF1 | 644 | KLF1 | 1641 |
| ctrl+MAP2K6 | 378 | MAP2K6 | 741 |
| ctrl+MAP7D1 | 358 | MAP7D1 | 635 |
| ctrl+MAPK1 | 381 | MAPK1 | 660 |
| ctrl+MEIS1 | 390 | MEIS1 | 682 |
| ctrl+OSR2 | 451 | OSR2 | 858 |
| ctrl+PRTG | 296 | PRTG | 581 |
| ctrl+PTPN12 | 164 | PTPN12 | 358 |
| ctrl+PTPN9 | 199 | PTPN9 | 363 |
| ctrl+RUNX1T1 | 356 | RUNX1T1 | 631 |
| ctrl+SAMD1 | 178 | SAMD1 | 430 |
| ctrl+SET | 278 | SET | 790 |
| ctrl+SLC38A2 | 280 | SLC38A2 | 280 |
| ctrl+SNAI1 | 160 | SNAI1 | 380 |
| ctrl+SPI1 | 336 | SPI1 | 506 |
| ctrl+STIL | 141 | STIL | 268 |
| ctrl+TBX2 | 295 | TBX2 | 564 |
| ctrl+TBX3 | 270 | TBX3 | 559 |
| ctrl+TGFBR2 | 328 | TGFBR2 | 537 |
| ctrl+UBASH3A | 315 | UBASH3A | 686 |
| ctrl+UBASH3B | 513 | UBASH3B | 983 |
| ctrl+ZBTB1 | 287 | ZBTB1 | 602 |
| ctrl+ZBTB25 | 247 | ZBTB25 | 590 |
