from pathlib import Path
import json,pandas as pd
root=Path(__file__).resolve().parents[1];p=root/'reports'
r=json.loads((p/'audit_evidence.json').read_text());q=json.loads((p/'perturbation_evidence.json').read_text());v=json.loads((p/'validation_evidence.json').read_text())
c=pd.read_csv(p/'perturbation_counts.csv',keep_default_na=False);can=pd.read_csv(p/'canonical_perturbation_counts.csv');targets=pd.read_csv(p/'target_genes.csv');bs=pd.read_csv(p/'barcode_suffix_counts.csv')
def table(df):
 cols=list(df.columns)
 return '\n'.join(['| '+' | '.join(cols)+' |','| '+' | '.join(['---']*len(cols))+' |']+['| '+' | '.join(str(x).replace('|','\\|') for x in row)+' |' for row in df.itertuples(index=False,name=None)])
s=f'''# AnnData data audit — Norman archive

Audited: 2026-09-11. **Completed; no model was trained.** The downloaded H5AD was read without changing it. Counts below describe this processed file, not necessarily the complete original experiment.

## Principal findings and evidence

| Requested finding | Conclusion | Direct evidence |
| --- | --- | --- |
| Physical observation | One cell-barcode-associated single-cell expression profile; interpreted as one captured cell. Singlet status is not independently established. | `obs.index.name == 'cell_barcode'`; 91,205 unique barcodes; one expression row and one condition assignment per barcode. Example: `AAACCTGAGGCATGTG-1`. |
| Variable | One measured gene-expression feature, identified by a gene ID and gene name. | `var.index.name == 'gene_id'`; 5,045 unique IDs and 5,045 unique `gene_name` values. Example: `ENSG00000187608` → `ISG15`. |
| Shape of X | **91,205 × 5,045**, cells × genes. | AnnData shape and HDF5 `X` shape both agree. CSR sparse matrix, `float32`. |
| Expression scale | **Natural-log1p-transformed, per-cell-scaled counts**, not raw counts or merely unlogged normalized counts. | All 37,317,477 nonzero entries were compared against the stored count layer; maximum absolute error in `X = log1p(count × scale[cell])` is {v['max_abs_log1p_scaled_count_error']:.9g}. |
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
- Archive members, read from its directory: `norman/`, `norman/go.csv`, and `norman/perturb_processed.h5ad`. The H5AD member was extracted to the requested path. The archive's H5AD-member CRC32 is **{2796850581:08x}**; reading the entire member through Python `zipfile` completed without a CRC error. Its uncompressed size is **2,228,610,012 bytes**. `go.csv` was not used to infer perturbation targets.
- Extracted H5AD SHA-256: `{r['sha256']}`. This records local file identity; no separately published SHA-256 was available for comparison.
- HDF5 top-level keys are exactly `X`, `layers`, `obs`, `uns`, `var`; the only stored layer is `counts`. There is no `raw`, `obsm`, `varm`, `obsp`, or processing-history object.
- Direct numerical scans covered the complete matrix and count layer. No fitting, prediction, dimensional reduction, differential-expression recomputation, or model training was performed.

## 1–2. Observations and variables

Each observation is a cell-barcode record with a transcript-expression vector and assigned perturbation. The index contains **91,205 unique, nonmissing barcodes**, each 18 characters long, including a numeric suffix. This is evidence for a single-cell profile as the unit of observation; it is not a perturbation average, a gene, or a biological replicate. The file provides no doublet score or singlet-validation field, so the physical claim is one *putative captured cell per barcode*, not a guarantee that every barcode arose from exactly one physical cell. Cell-line identity remains unresolved as noted above.

Each variable is an annotated gene feature. IDs have the `ENSG` naming pattern (for example `ENSG00000187608`), paired with gene names (for example `ISG15`). These are gene-expression features, not guide sequences: `var` contains only `gene_name`, and the matrix associates their expression with cell barcodes. Gene names include protein-coding-style symbols and noncoding-style names such as `LINC01342` and `RP11-34P13.8`; no biotype annotation is supplied. Both feature identifiers and names are unique. The file does not record the selection rule that produced these 5,045 features; it cannot be assumed that they are all genes originally measured or all highly variable genes.

Evidence: [complete observation annotations](obs_metadata.csv), [complete feature annotations](var_metadata.csv), [structural and numerical evidence](audit_evidence.json).

## 3–4. Matrix shape, counts, and transformation

`X` has **460,129,225 possible entries**, of which **37,317,477 are nonzero** ({100*37317477/460129225:.4f}% density; {100*(1-37317477/460129225):.4f}% zeros). It has no empty rows or empty columns. All values are finite and nonnegative. Implicit sparse zeros are included in these conclusions.

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

All **37,317,477** nonzero entries match within **2 × 10⁻⁶** absolute error; the worst error is **{v['max_abs_log1p_scaled_count_error']:.9g}**. Zero locations match exactly. The per-cell scales range from **0.2474850895 to 6.5645891100**, with median **1.0000000355**.

For the first barcode, `AAACCTGAGGCATGTG-1`, directly observed entries are:

| Gene | Stored count | X | `expm1(X)/count` |
| --- | --- | --- | --- |
| ISG15 | 1 | 0.6815017462 | 0.9768442225 |
| ENO1 | 4 | 1.5907396078 | 0.9768442802 |
| SRM | 2 | 1.0830547810 | 0.9768443278 |

Thus the data support **per-cell scaling followed by natural log1p**. They do not establish a particular normalization function, denominator gene set, or target total. The retained genes' `expm1(X)` row sums range from **1,964.0674 to 9,479.6816** (median **2,963.7360**), so those sums are not a constant 10,000. No `uns['log1p']` or normalization settings are stored. A tested 10,000-total hypothesis does not recover consistently integer original totals (only **{100*v['normalization_10000_hypothesis']['fraction_within_0_01_of_integer']:.3f}%** of `10000/scale` values lie within 0.01 of an integer); this does not reconstruct or prove the original normalization recipe.

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

{table(bs.rename(columns={'barcode_suffix':'Barcode suffix','cells':'Cells','control_cells':'Controls','condition_labels':'Observed condition labels'}))}

All **284 condition labels occur in every suffix group**, and controls occur in every group. This is evidence of eight barcode partitions and is consistent with multiple source libraries or groups, but the file does not map suffixes to collection dates, lanes, libraries, or biological/technical replicates. Therefore **“no batches exist” is unsupported**, and **“eight experimentally verified batches” would also overstate the evidence**. The suffix is a candidate technical grouping, not an authenticated batch annotation. The audit did not test expression batch effects.

The complete barcodes are unique. Removing suffixes introduces **4,768 duplicate occurrences beyond the first**, so suffixes must be preserved as part of cell identity. Evidence: [suffix summary](barcode_suffix_counts.csv), [condition × suffix counts](condition_by_barcode_suffix.csv), [validation_evidence.json](validation_evidence.json).

## 10. Genes actually targeted according to the assignments

The **105 distinct symbols below** come exclusively from observed perturbation labels, after removing `ctrl`. All have at least one single-target condition in this file and all match `var['gene_name']` exactly. The other **4,940 measured features** are not named as targets in any observed condition. A measured gene or a gene in a stored ranking was not counted as a target on that basis.

'''
for i in range(0,len(q['target_genes']),15):s+=', '.join('`'+g+'`' for g in q['target_genes'][i:i+15])+'.\n\n'
s+='''The symbols are retained as stored, including names such as `C19orf26`, `KIAA1804`, and `RHOXF2BB`; no external alias conversion was guessed. [target_genes.csv](target_genes.csv) gives each target's exact feature ID, number of cells in single-target conditions, number in two-target conditions, total assigned cells, and number of original labels. Per-gene totals overlap for two-target cells and must not be summed as a count of unique cells.

“Targeted” here means **assigned an annotated perturbation against that gene**. Successful molecular activation, guide identity, and off-target effects cannot be established from the supplied fields.

## Reproducibility and evidence files

The audit can be rerun from the project root, with the H5AD already present:

```powershell
python scripts/inspect_norman.py
python scripts/audit_perturbations.py
python scripts/validate_norman.py
python scripts/write_data_audit.py
```

'''
s+='Inspected with '+', '.join(f'`{k} {ver}`' for k,ver in r['versions'].items())+'.\n\n'
s+='''The scripts inspect stored data and perform deterministic summaries and numerical identity checks; they do not train models. Numerical evidence is retained in [audit_evidence.json](audit_evidence.json), [perturbation_evidence.json](perturbation_evidence.json), and [validation_evidence.json](validation_evidence.json). Complete `obs` and `var` exports and all count tables are linked above. The download archive's member manifest is retained in `data/norman/archive_manifest.json`.

## Appendix A. Every original perturbation label and cell count

The target-set total is repeated when two original labels map to the same target set. **Sum the original-label cells column to obtain 91,205; do not sum the repeated target-set totals.** Use the separate canonical table for unique target-set totals.

'''
c['target_set_cells']=c.canonical_perturbation.map(can.set_index('canonical_perturbation').n_cells)
s+=table(c[['condition','n_cells','canonical_perturbation','target_set_cells']].rename(columns={'condition':'Original condition','n_cells':'Cells','canonical_perturbation':'Target set','target_set_cells':'Target-set cells'}))+'\n'
(p/'DATA_AUDIT.md').write_text(s,encoding='utf-8')
assert c.n_cells.sum()==91205 and can.n_cells.sum()==91205
assert len(c)==284 and len(can)==237 and len(targets)==105 and targets.measured_in_X.all()
assert sum(q['cell_counts_by_n_targets'].values())==91205
assert len(c.condition.unique())==284
print('Wrote',p/'DATA_AUDIT.md',len(s),'characters; all count-total checks passed.')
