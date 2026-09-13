"""Tiny tests for the ranking-plot layout helper and report splicer."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from finalize_modeling_report import splice_report
from modeling_analysis import plot_pair_error_distribution_ranking, _try_pyplot


class TestPairRankingLayout(unittest.TestCase):
    def test_constrained_layout_keeps_two_panels(self):
        plt = _try_pyplot()
        if plt is None:
            self.skipTest("matplotlib not importable")
        frame = pd.DataFrame(
            {
                "canonical_perturbation": [f"G{i}+H{i}" for i in range(8)],
                "additive_rmse": [0.01 * (i + 1) for i in range(8)],
            }
        )
        fig = plot_pair_error_distribution_ranking(frame, plt)
        self.assertTrue(bool(fig.get_constrained_layout()))
        self.assertEqual(len(fig.axes), 2)
        plt.close(fig)


class TestReportSplicer(unittest.TestCase):
    def test_replaces_existing_marker_block(self):
        original = "before\n<!-- BEGIN_AUDITED_SUMMARY -->\nold\n<!-- END_AUDITED_SUMMARY -->\nafter\n"
        out = splice_report(original, "## Audited\nnew")
        self.assertIn("new", out)
        self.assertNotIn("old", out)
        self.assertEqual(out.count("BEGIN_AUDITED_SUMMARY"), 1)


if __name__ == "__main__":
    unittest.main()
