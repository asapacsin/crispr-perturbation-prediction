# Worker resource review (not a final scientific report)

Orchestrator follow-up: directly streamed the archive member without extracting it. `norman/go.csv` has 673,899 source/target/importance edges and 3,127 nodes, all measured genes, but no named term memberships. See modeling_resource_evidence.json. Lane-balanced aggregation is compatible with split B; globally holding out each pair is intentional. Equal lane weighting is a declared new estimand, not evidence that the old split is invalid. The new modeling_protocol.yaml sets MSE as the only primary metric; historical proposal fields are not the active protocol.

I'll do a bounded read-only pass: the provenance/methods/schema files first, then local GO/gene-set resources and existing H5AD/plotting helpers—no training or dense-matrix loads.`METHODS_DECISION_BRIEF.md` isn't at the repo root. Next I’ll locate that brief, `analyze_provenance.py`, local GO files, and H5AD/plotting helpers.`go.csv` is documented as a zip member, not a loose file. I’ll inspect the archive listing, GEARS GO loaders, H5AD helpers, and lane/interaction notes without extracting or loading dense `X`.Next I’ll peek `go.csv` inside the zip (header only), check GEARS graph construction, and inspect installed linear/plotting packages plus lane-composition tables.`METHODS_DECISION_BRIEF.md` is at `reports/METHODS_DECISION_BRIEF.md`, not the repo root. Shell was blocked, so `norman/go.csv` inside the zip was not streamed and current `import` status for plotting/linear packages was not rechecked.

## 1. Contradictions / pitfalls for lane-balanced single-effect and interaction analysis

**Lanes are technical, not biological replicates.** Schema and provenance treat the eight barcode suffixes as 10x Chromium lanes of one pooled K562 CRISPRa experiment, not donors. Using \(n=8\) as biological replication, or treating protocol D as a donor holdout, is the documented error.

```312:315:config/dataset_schema.yaml
groups:
  n: 8
  barcode_suffix_matches_source_gemgroup: true
  interpretation: technical_10x_chromium_lanes_not_biological_replicates_or_donors
```

`analyze_provenance.py` already finds a **condition × lane association** (χ² \(p=1.6\times10^{-4}\), Cramér’s \(V\approx0.059\)) and small but significant **control-only** Kruskal lane differences (ε²≈0.004). Control mean profiles still correlate 0.9983–0.9998. Every lane has all 284 labels and all 105 targets, but **counts are not exchangeable**.

```19:36:scripts/analyze_provenance.py
ct=pd.crosstab(o.condition,o.barcode_group);chi,pval,df,expected=st.chi2_contingency(ct)
# ...
# Control-only summaries avoid perturbation-composition confounding; no causal batch-effect claim.
```

Per-lane cell totals are 10,606–12,409; controls 821–1,008 (`reports/barcode_group_summary.csv`). Some **gene × lane** cells are too thin for a stable per-lane single-effect (e.g. ARID1A in lane 1: 24 single cells, 0 doubles in `reports/target_gene_by_barcode_group.csv`). Condition sizes also span 49–1,641 noncontrol cells. A pooled single-effect is **cell-weighted**; a lane-equal mean is a different estimand. Systema’s cell-weighted perturbed mean is **not** the brief’s condition-balanced centroid.

**Protocol B is not a lane-balanced interaction design.** Pair assignments are global hashes; **all lanes of a held-out pair go to the same split**. Only controls are split 10/10/rest **inside each lane**.

```55:60:scripts/analyze_provenance.py
pairs=rank(sorted(o.loc[o.n_targets==2,'canonical_perturbation'].unique()),'B_pairs');assert len(pairs)==131
pairmap={p0:('test' if k<26 else 'validation' if k<52 else 'train') for k,p0 in enumerate(pairs)}
manifest['B_heldout_combination']=o.canonical_perturbation.map(pairmap).fillna('train')
```

B tests **predictive combination generalization**, not a causal GI or paired treatment effect (`reports/METHODS_DECISION_BRIEF.md` §5). Additive \(\delta_i+\delta_j\) on stored **mean-log1p** is an assumption, not a count-space GI. If singles and pairs have different lane weights, \(\mu_{ij}-(\delta_i+\delta_j)\) can be **lane-mix residual**, not interaction. Ahlmann-Eltze (cited in the brief): on Norman doubles, interaction-vs-additive was no better than no-change.

**Document contradictions that affect this analysis**

| Conflict | What it implies |
| --- | --- |
| Schema lists **two** primaries (`macro_condition_mse`, `macro_condition_mae`); brief insists on **one** winner (MSE), MAE secondary | `config/dataset_schema.yaml` 430–432 vs `reports/METHODS_DECISION_BRIEF.md` 96–101 |
| Schema skill pair is `common_shift` + `additive`; brief skill pair is **condition-balanced perturbed mean** + additive (no-change is not the justifying pair) | same lines |
| Historical `reports/DATA_AUDIT.md` 17: “experimental batch identities are not annotated” | Later provenance **did** resolve suffixes as GEM lanes. Do not treat the audit sentence as current design. |
| `orientation_equivalent_for_target_set: true` vs `experimental_orientation_equivalence: unknown` | Collapsing `GENE+ctrl` / `ctrl+GENE` (e.g. KLF1 997+644) is a **split convenience**, not an experimental result (`config/dataset_schema.yaml` 75–76; `reports/DATA_AUDIT.md` 119) |
| Proposed delta uses **training-split controls only**, not per-lane controls | Mixing per-lane control centering with the schema formula is a different estimand (`config/dataset_schema.yaml` 420–424) |
| Same-reference MSE identity | \(\mathrm{MSE}(\mu-r,\hat\mu-r)=\mathrm{MSE}(\mu,\hat\mu)\) for shared \(r\); “delta MSE” is not a new number (`reports/METHODS_DECISION_BRIEF.md` 82) |
| Target-excluded 4,940 drops **all 105** panel genes | Removes other CRISPRa program TFs, not only the queried loci (brief §3) |
| Unpaired design | No cell-level counterfactual GI (`pairing` in schema 375–380) |
| GEARS README: singles-only training **cannot** reliably predict combos | `reports/provenance_sources/GEARS/README.md` 23; `filter_pert_in_go` **drops** conditions whose tokens are off the GO graph (`utils.py` 318–337) |

HVG was selected on the full object; stored `uns` DE is banned for fitting. `bulk_fitness` is leakage if used as a covariate.

---

## 2. Local gene-set / GO resources

**`data/norman/go.csv` is not extracted.** File-tool read of that path failed. `scripts/download_norman.py` 29–31 extracts **only** `perturb_processed.h5ad`.

The file exists as a **zip member**:

```8:13:data/norman/archive_manifest.json
    "name": "norman/go.csv",
    "bytes": 19700177,
    "compressed_bytes": 3866939,
    "crc": 3995457773
```

`reports/DATA_AUDIT.md` 27: archive members are `norman/`, `norman/go.csv`, `norman/perturb_processed.h5ad`; **`go.csv` was not used to infer targets**.

**Expected schema (GEARS writer), not a counted local table** — zip body was not opened (shell blocked):

```126:154:reports/provenance_sources/GEARS/gears/utils.py
    go_path = os.path.join(data_path, data_name, 'go.csv')
    if os.path.exists(go_path):
        return pd.read_csv(go_path)
    # ...
        df_edge_list = pd.DataFrame(further_filter).rename(columns = {0: 'gene1',
                                                                      1: 'gene2',
                                                                      2: 'score'})
        df_edge_list = df_edge_list.rename(columns = {'gene1': 'source',
                                                      'gene2': 'target',
                                                      'score': 'importance'})
        df_edge_list.to_csv(go_path, index = False)
```

Columns: **`source`, `target`, `importance`**. `importance` is **Jaccard overlap of GO-term ID sets** between two genes (keep score > 0.1). There are **no GO term IDs, names, or definitions** in this CSV.

**Coverage is a perturbation-graph, not a 5,045-gene response GMT.**

- `get_go_auto` builds edges only among `gene_list` keys present in `gene2go`.
- `make_GO` subsets `gene2go` to `pert_list` (`utils.py` 207–218).
- Runtime default GO graph is **not** this zip file: `get_similarity_network(..., network_type='go')` downloads `go_essential_all/go_essential_all.csv` (Dataverse 6934319) unless `default_pert_graph=False` (`utils.py` 247–257).
- `PertData.set_pert_genes` default loads `essential_all_data_pert_genes.pkl`; only the non-default branch uses `adata.var['gene_name']` plus pert tokens (`pertdata.py` 104–128). Nodes missing from `gene2go` are dropped; GEARS then **filters those conditions out of the AnnData** (`pertdata.py` 184–194).
- `gene2go.pkl` / `gene2go_all.pkl`, GMT/OBO/GAF, and GEARS `genes_with_hi_mean.npy` are **absent** locally (glob found none).
- Stored aliases (`RHOXF2BB`, `C19orf26`, `C3orf72`, `KIAA1804`) are **not** HGNC-equivalent (`config/dataset_schema.yaml` 81–83). They are the tokens most likely to miss standard GO keys.
- Measured response genes are 5,045 `var.gene_name` values (`reports/var_metadata.csv`: `gene_id,gene_name`). All **105 assigned targets** are among them (`reports/DATA_AUDIT.md` 18; `reports/target_genes.csv` `measured_in_X`). That does **not** imply `go.csv` covers the other ~4,940 genes.

**GO term descriptions: no.** Local GO artifacts are Jaccard **edges**, not a term thesaurus. `gene2go_*` would be gene → term-ID sets only, and those pickles are not in the tree.

GEARS GI scoring can restrict to a high-mean gene list (`get_GI_genes_idx`, `utils.py` 534–548); that `.npy` is not in this repo. `GI_genes_file=None` uses all measured genes (`gears.py` 405–409).

---

## 3. H5AD helpers and plotting / linear stack

**Prefer these if the goal is annotations without dense `counts`:**

| Helper | What it loads | Path |
| --- | --- | --- |
| `load_obs_frame` / `load_var_frame` / `matrix_shape` | HDF5 `obs`/`var` + shapes only; **does not load counts payload** | `scripts/validate_dataset_schema.py` 129–183, 220–228 |
| `canonicalize_condition`, split checks | labels only | `scripts/dataset_schema_lib.py` |
| `audit_perturbations.py` | `reports/obs_metadata.csv`, `reports/var_metadata.csv` | no H5AD |

**Full AnnData loads (includes `X` / can touch `counts`):** `scripts/inspect_norman.py` 4; `scripts/validate_norman.py` 3; `scripts/analyze_provenance.py` 4 (`ad.read_h5ad` then `a.layers['counts'].sum` and `a.X` control means). GEARS `PertData.load` uses `sc.read_h5ad` (`pertdata.py` 170–171).

**Recorded project stack** (`reports/audit_evidence.json` 5–11; `AGENTS.md` 15): `anndata 0.13.3.post0`, `numpy 2.5.2`, `scipy 1.18.1`, `pandas 3.0.5`, `h5py 3.16.0`. Schema validator **requires PyYAML** (`validate_dataset_schema.py` 17–23). Linear summaries already used: `scipy.stats` χ²/Kruskal (`analyze_provenance.py` 2, 19–35).

**Not in the project’s recorded inspect list:** matplotlib, seaborn, scikit-learn, statsmodels, scanpy, plotly. Those appear only in **reference** GEARS (`requirements.txt`: sklearn, scanpy, scipy; `gears.py` 431–432 imports seaborn/matplotlib; `utils.py` 12 `TheilSenRegressor`). No root `requirements.txt`. **Could not re-`import` the live env** (shell blocked).

**Could not inspect:** unique genes / overlap of zip `norman/go.csv` vs 5,045 `gene_name`; whether `gene2go` would include `RHOXF2BB`; current presence of matplotlib/sklearn in the interpreter.