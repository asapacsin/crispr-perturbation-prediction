"""Build reports/verified_target_guide_sequences.csv from local source tables.

Does not guess missing names. Does not copy bulk_fitness into the verified table.
Compares protospacers for every gene shared by both author tables and requires
a 20-base ACGT sequence wherever a sequence is present.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
REPORTS = ROOT / "reports"
SOURCES = REPORTS / "provenance_sources"
OUT = REPORTS / "verified_target_guide_sequences.csv"

RHOXF2_SEQUENCE = "GCATGCGTTGTCCTGTAGAC"
RHOXF2_ALIASES = "RHOXF2;RHOXF2B;RHOXF2BB"
ACGT20 = re.compile(r"^[ACGT]{20}$")


def assert_acgt20(sequence: str, gene: str, table: str) -> None:
    if not ACGT20.fullmatch(sequence):
        raise SystemExit(
            f"{table} protospacer for {gene} is not 20 ACGT characters: {sequence!r}"
        )


def compare_shared_sequences(sharing: pd.DataFrame, emaps: pd.DataFrame) -> None:
    sharing_by_gene = sharing.set_index("gene", verify_integrity=True)
    emaps_by_gene = emaps.set_index("gene", verify_integrity=True)
    mismatches = []
    for gene in sorted(set(sharing_by_gene.index) & set(emaps_by_gene.index)):
        share_seq = str(sharing_by_gene.loc[gene, "protospacer"]).upper()
        emap_seq = str(emaps_by_gene.loc[gene, "protospacer"]).upper()
        assert_acgt20(share_seq, gene, "data_sharing")
        assert_acgt20(emap_seq, gene, "final_emaps")
        if share_seq != emap_seq:
            mismatches.append(gene)
    sharing_only = set(sharing_by_gene.index) - set(emaps_by_gene.index)
    for gene in sorted(sharing_only):
        assert_acgt20(str(sharing_by_gene.loc[gene, "protospacer"]).upper(), gene, "data_sharing")
    emaps_only = set(emaps_by_gene.index) - set(sharing_by_gene.index)
    for gene in sorted(emaps_only):
        assert_acgt20(str(emaps_by_gene.loc[gene, "protospacer"]).upper(), gene, "final_emaps")
    if mismatches:
        raise SystemExit(
            "protospacer mismatch on genes present in both guide tables: "
            + ",".join(mismatches)
        )


def main() -> Path:
    targets = pd.read_csv(REPORTS / "target_genes.csv")
    sharing = pd.read_csv(SOURCES / "guide_sequences_data_sharing.csv")
    emaps = pd.read_csv(SOURCES / "guide_sequences_final_emaps.csv")
    compare_shared_sequences(sharing, emaps)
    sharing_by_gene = sharing.set_index("gene", verify_integrity=True)
    emaps_by_gene = emaps.set_index("gene", verify_integrity=True)

    rows = []
    for row in targets.itertuples(index=False):
        stored = row.target_gene
        lookup = stored
        mapping_rule = "exact_source_gene_name"
        source_aliases = stored
        unresolved = "resolved"
        notes = "Exact gene-name match in at least one author guide table."

        if stored == "RHOXF2BB":
            lookup = "RHOXF2"
            mapping_rule = "documented_source_RHOXF2_to_stored_RHOXF2BB"
            source_aliases = RHOXF2_ALIASES
            unresolved = "resolved_dataset_local_only"
            notes = (
                "Source guide_identity uses RHOXF2. Author GI_generate_populations "
                "cell 3 relabels RHOXF2 to RHOXF2B because the guide lands between "
                "RHOXF2/RHOXF2B with a stronger reported effect on B. GEARS Norman19 "
                "notebook maps var name RHOXF2B to stored RHOXF2BB. data_sharing lists "
                "RHOXF2 and RHOXF2B with the same protospacer; final_emaps lists only "
                "RHOXF2. This is dataset-local label reconciliation, not gene-symbol "
                "equivalence or molecular-specificity proof."
            )
        elif stored not in sharing_by_gene.index and stored not in emaps_by_gene.index:
            lookup = ""
            mapping_rule = "none"
            source_aliases = stored
            unresolved = "unresolved_no_source_sequence"
            notes = (
                "Stored target is measured and assigned in the H5AD, but neither "
                "author guide-sequence table contains this gene name. No sequence "
                "was guessed."
            )

        sharing_hit = lookup in sharing_by_gene.index if lookup else False
        emaps_hit = lookup in emaps_by_gene.index if lookup else False
        extra_alias_row = stored == "RHOXF2BB" and "RHOXF2B" in sharing_by_gene.index
        guide_id = ""
        sequence = ""
        if sharing_hit:
            guide_id = str(sharing_by_gene.loc[lookup, "guide_id"])
            sequence = str(sharing_by_gene.loc[lookup, "protospacer"]).upper()
        elif emaps_hit:
            guide_id = str(emaps_by_gene.loc[lookup, "guide_id"])
            sequence = str(emaps_by_gene.loc[lookup, "protospacer"]).upper()
        if sequence:
            assert_acgt20(sequence, stored, "verified_mapping")

        if extra_alias_row:
            alias_seq = str(sharing_by_gene.loc["RHOXF2B", "protospacer"]).upper()
            assert_acgt20(alias_seq, "RHOXF2B", "data_sharing")
            if sequence and alias_seq != sequence:
                raise SystemExit(
                    "RHOXF2 and RHOXF2B protospacers differ; refuse to collapse them"
                )
            notes += " Both RHOXF2 and RHOXF2B data_sharing rows share protospacer GCATGCGTTGTCCTGTAGAC."

        rows.append(
            {
                "stored_target_gene": stored,
                "stored_gene_id": row.gene_id,
                "source_gene": lookup,
                "source_aliases": source_aliases,
                "source_guide_id": guide_id,
                "protospacer": sequence,
                "present_in_data_sharing": sharing_hit or extra_alias_row,
                "present_in_final_emaps": emaps_hit,
                "mapping_rule": mapping_rule,
                "unresolved_status": unresolved,
                "notes": notes,
            }
        )

    out = pd.DataFrame(rows)
    if len(out) != 105:
        raise SystemExit(f"expected 105 stored targets, wrote {len(out)}")
    if (out["stored_target_gene"] == "SLC38A2").sum() != 1:
        raise SystemExit("SLC38A2 missing from verified table")
    slc = out.loc[out["stored_target_gene"] == "SLC38A2"].iloc[0]
    if slc["protospacer"] or slc["unresolved_status"] != "unresolved_no_source_sequence":
        raise SystemExit("SLC38A2 must remain unresolved without a guessed sequence")
    rho = out.loc[out["stored_target_gene"] == "RHOXF2BB"].iloc[0]
    if rho["source_gene"] != "RHOXF2" or rho["protospacer"] != RHOXF2_SEQUENCE:
        raise SystemExit("RHOXF2BB mapping failed")
    if "bulk_fitness" in out.columns:
        raise SystemExit("bulk_fitness must not appear in the verified mapping table")
    out.to_csv(OUT, index=False)
    print(f"wrote {OUT} rows={len(out)}")
    return OUT


if __name__ == "__main__":
    main()
