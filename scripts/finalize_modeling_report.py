"""Write the audited biological summary from frozen diagnostic tables.

Does not refit, reload X, or recompute long effect tables.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

from modeling_analysis import refresh_pair_error_distribution_ranking

DIAG = ROOT / "reports" / "modeling" / "diagnostics"
REPORT = ROOT / "reports" / "MODELING_REPORT.md"
SUMMARY = ROOT / "reports" / "modeling" / "AUDITED_SUMMARY.md"
BEGIN = "<!-- BEGIN_AUDITED_SUMMARY -->"
END = "<!-- END_AUDITED_SUMMARY -->"

# Destination-relative hrefs. AUDITED_SUMMARY lives in reports/modeling/;
# MODELING_REPORT lives in reports/.
LINKS = {
    "summary": {
        "diagnostics": "../MODELING_DIAGNOSTICS.md",
        "root_verification": "root_verification.json",
        "ranking": "describe/figures/pair_error_distribution_ranking.png",
        "heatmap": "describe/figures/strongest_interaction_heatmap.png",
    },
    "report": {
        "diagnostics": "MODELING_DIAGNOSTICS.md",
        "root_verification": "modeling/root_verification.json",
        "ranking": "modeling/describe/figures/pair_error_distribution_ranking.png",
        "heatmap": "modeling/describe/figures/strongest_interaction_heatmap.png",
    },
}


def _fmt(value, digits=6) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return "NA"
    if isinstance(value, float):
        return f"{value:.{digits}g}"
    return str(value)


def _md_table(rows: list[dict], columns: list[str]) -> str:
    header = "| " + " | ".join(columns) + " |"
    sep = "| " + " | ".join("---" for _ in columns) + " |"
    body = [
        "| " + " | ".join(_fmt(row.get(col, "")) for col in columns) + " |"
        for row in rows
    ]
    return "\n".join([header, sep, *body])


def build_audited_summary(destination: str = "summary") -> str:
    genes = pd.read_csv(DIAG / "example_pair_genes.csv")
    own = json.loads((DIAG / "own_target_counts.json").read_text(encoding="utf-8"))
    shift = json.loads((DIAG / "shared_shift_summary.json").read_text(encoding="utf-8"))
    sizes = pd.read_csv(ROOT / "reports" / "modeling" / "prepare" / "condition_sizes.csv")
    pairs = pd.read_csv(ROOT / "reports" / "modeling" / "describe" / "pair_summaries.csv")
    n_pair = int(sizes.loc[sizes["canonical_perturbation"] == "CEBPA+CEBPB", "total_n"].iloc[0])
    cebpa = genes.loc[genes["canonical_perturbation"] == "CEBPA+CEBPB"].head(3)
    zc3 = pairs.loc[pairs["canonical_perturbation"] == "CEBPA+ZC3HAV1"].iloc[0]
    lst1 = cebpa.iloc[0]
    table = _md_table(
        [
            {
                "gene": row.gene_name,
                "A": row.effect_a,
                "B": row.effect_b,
                "AB": row.observed_ab,
                "A+B": row.additive,
                "residual": row.residual,
                "SE": row.conditional_cell_sampling_se,
            }
            for row in cebpa.itertuples(index=False)
        ],
        ["gene", "A", "B", "AB", "A+B", "residual", "SE"],
    )
    links = LINKS[destination]
    return f"""## Audited biological summary

Numbers below are copied from [MODELING_DIAGNOSTICS.md]({links["diagnostics"]}) and `reports/modeling/diagnostics/`. They did not enter fit or alpha selection.

**Plain-language example.** On stored log1p `X`, CEBPA alone raises LST1 by {_fmt(lst1.effect_a)}, CEBPB alone by {_fmt(lst1.effect_b)}, and the pair reaches {_fmt(lst1.observed_ab)} rather than the additive {_fmt(lst1.additive)}. The LST1 residual is {_fmt(lst1.residual)} (conditional SE {_fmt(lst1.conditional_cell_sampling_se)}). Root recomputed that same LST1 contrast on original stored log1p `X` ([root_verification.json]({links["root_verification"]})). This is one-experiment cell-sampling uncertainty, not a biological p-value. The pair has **n={n_pair}** cells. Construct-orientation alternative additive RMSE for this pair is **0.09046–0.12455**, so a residual remains but its magnitude is orientation-sensitive. Do not read this as proven synergy.

Top 3 |residual| genes in `CEBPA+CEBPB`:

{table}

**Own-target check.** {own["n_positive"]} of {own["n_targets"]} measured targets have a positive own-locus log1p effect ({own["n_negative"]} negative). That is a normalized consistency check, not molecular validation of guides or protein.

**A more additive contrast.** `CEBPA+ZC3HAV1` has residual/response RMS **{_fmt(zc3.residual_over_response_rms, 3)}** (observed RMS {_fmt(zc3.observed_rms)}, residual RMS {_fmt(zc3.additive_rmse)}). That is closer to additivity than CEBPA+CEBPB and is **not** perfect additivity.

**Modules, not new GO.** Named GO enrichment is unavailable (graph has source/target/importance only). Data-driven SVD PC1 separates LST1/AIF1/CSF3R from HBZ/HBG1/HBG2. Norman et al. 2019 already interpret related programs as granulocyte- versus erythroid-like; that is a **supported existing interpretation**, not a new enrichment test.

**Shared shift.** Median residual squared-norm fraction along the mean training single-delta is **{_fmt(shift["median_squared_norm_fraction"], 4)}**. Alignment does not exclude other confounding.

**Held-out prediction (already frozen).** Ridge α=0.1 test macro MSE **0.001599794** vs additive **0.002372021** (32.56% lower; 15/26 pairs). The 26 test pairs are consumed. Next models need nested resampling or new external holdout.

Full diagnostic report: [MODELING_DIAGNOSTICS.md]({links["diagnostics"]}). Figures: [pair ranking]({links["ranking"]}), [zero-centered interaction heatmap]({links["heatmap"]}).
"""


def splice_report(report_text: str, section: str) -> str:
    block = f"{BEGIN}\n{section.rstrip()}\n{END}\n"
    if BEGIN in report_text and END in report_text:
        start = report_text.index(BEGIN)
        stop = report_text.index(END) + len(END)
        return report_text[:start] + block + report_text[stop:]
    anchor = "## Repeatability"
    if anchor in report_text:
        return report_text.replace(anchor, block + "\n" + anchor, 1)
    return report_text.rstrip() + "\n\n" + block


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--refresh-ranking-plot",
        action="store_true",
        help="Also regenerate pair_error_distribution_ranking.png from pair_summaries.csv",
    )
    args = parser.parse_args()
    if not (DIAG / "example_pair_genes.csv").is_file():
        raise SystemExit("diagnostics tables are not ready")
    summary_section = build_audited_summary("summary")
    SUMMARY.write_text(summary_section + "\n", encoding="utf-8")
    if REPORT.is_file():
        REPORT.write_text(
            splice_report(REPORT.read_text(encoding="utf-8"), build_audited_summary("report")),
            encoding="utf-8",
        )
    if args.refresh_ranking_plot:
        path = refresh_pair_error_distribution_ranking()
        print("finalize: refreshed " + str(path), flush=True)
    print("finalize: wrote " + str(SUMMARY), flush=True)
    print("finalize: updated " + str(REPORT), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
