from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from spacer.html_report import write_html_report


class HtmlReportTests(unittest.TestCase):
    def test_writes_self_contained_interactive_report(self) -> None:
        analysis = [
            {
                "summary_level": "run", "run_basename": "run_1", "run_count": "1",
                "total_ms2": "10", "scorable_ms2": "8", "low_interference": "5",
                "likely_chimeric": "3", "indeterminate": "2",
                "likely_chimeric_fraction": "0.375", "indeterminate_fraction": "0.2",
                "replicate_likely_fraction_median": "", "replicate_likely_fraction_min": "",
                "replicate_likely_fraction_max": "", "interference_threshold": "0.3",
            },
            {
                "summary_level": "single_group", "run_basename": "single_group", "run_count": "1",
                "total_ms2": "10", "scorable_ms2": "8", "low_interference": "5",
                "likely_chimeric": "3", "indeterminate": "2",
                "likely_chimeric_fraction": "0.375", "indeterminate_fraction": "0.2",
                "replicate_likely_fraction_median": "0.375", "replicate_likely_fraction_min": "0.375",
                "replicate_likely_fraction_max": "0.375", "interference_threshold": "0.3",
            },
        ]
        agreement = [
            {
                "run_basename": "run_1", "subset": "all_matched", "score_rows": "10",
                "matched_finite": "7", "pd_interference_missing": "1",
                "spacer_indeterminate": "2", "scoring_without_pd_spectrum": "0",
                "pearson_r": "0.8", "spearman_rho": "0.7",
                "median_signed_difference_pp": "-2", "median_absolute_difference_pp": "4",
            }
        ]
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "report.html"
            write_html_report(
                path,
                analysis_rows=analysis,
                sensitivity_rows=[],
                interference_rows=[
                    {
                        "summary_level": "single_group",
                        "run_basename": "single_group",
                        "scorable_ms2": "8",
                        "interference_band": "0_to_under_10_percent",
                        "lower_inclusive_percent": "0",
                        "upper_inclusive_percent": "10",
                        "scan_count": "2",
                        "fraction_of_scorable_ms2": "0.25",
                        "median_interference_fraction": "0.2",
                        "p10_interference_fraction": "0.01",
                        "p25_interference_fraction": "0.1",
                        "p75_interference_fraction": "0.4",
                        "p90_interference_fraction": "0.8",
                    }
                ],
                agreement_rows=agreement,
                discordant_rows=[],
                validation_rows=[
                    {
                        "run_basename": "run_1", "component": "score_rows",
                        "severity": "info", "status": "pass", "detail": "ok",
                    }
                ],
                plot_rows=[],
            )
            page = path.read_text(encoding="utf-8")
        self.assertIn("Spacer Single-group report", page)
        self.assertIn("PD values are descriptive reference context only", page)
        self.assertIn("const data =", page)
        self.assertIn("Discordant scans", page)
        self.assertIn("Interference degree", page)
        self.assertIn("PD-matched plots", page)
