from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from spacer.reporting import ReportingError, build_analysis, summarize_interference_distribution


def _write(path: Path, header: str, rows: list[str]) -> None:
    path.write_text(header + "\n" + "\n".join(rows) + "\n", encoding="utf-8")


class ReportingTests(unittest.TestCase):
    def test_builds_descriptive_group_and_sensitivity_summaries(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            _write(
                directory / "run_summary.tsv",
                "\t".join(
                    [
                        "run_basename", "total_ms2", "scorable_ms2", "low_interference",
                        "likely_chimeric", "indeterminate", "likely_chimeric_fraction",
                        "interference_threshold",
                    ]
                ),
                [
                    "run_1\t10\t8\t5\t3\t2\t0.375\t0.3",
                    "run_2\t20\t18\t9\t9\t2\t0.5\t0.3",
                ],
            )
            _write(
                directory / "sensitivity_summary.tsv",
                "\t".join(
                    [
                        "run_basename", "interference_threshold", "total_ms2", "scorable_ms2",
                        "low_interference", "likely_chimeric", "indeterminate",
                        "likely_chimeric_fraction", "minimum_competitor_prior_ms1_detections",
                    ]
                ),
                [
                    "run_1\t0.2\t10\t8\t4\t4\t2\t0.5\t2",
                    "run_2\t0.2\t20\t18\t8\t10\t2\t0.555555555556\t2",
                    "run_1\t0.3\t10\t8\t5\t3\t2\t0.375\t2",
                    "run_2\t0.3\t20\t18\t9\t9\t2\t0.5\t2",
                ],
            )
            analysis, sensitivity, report = build_analysis(directory)
            self.assertEqual(3, len(analysis))
            self.assertEqual("single_group", analysis[-1]["summary_level"])
            self.assertEqual("0.461538461538", analysis[-1]["likely_chimeric_fraction"])
            self.assertEqual(["0.2", "0.3"], [row["interference_threshold"] for row in sensitivity])
            self.assertIn("no between-method comparison", report)

    def test_rejects_inconsistent_score_arithmetic(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            header = "\t".join(
                [
                    "run_basename", "total_ms2", "scorable_ms2", "low_interference",
                    "likely_chimeric", "indeterminate", "likely_chimeric_fraction",
                    "interference_threshold",
                ]
            )
            _write(directory / "run_summary.tsv", header, ["run_1\t10\t8\t5\t2\t2\t0.25\t0.3"])
            _write(
                directory / "sensitivity_summary.tsv",
                header + "\tminimum_competitor_prior_ms1_detections",
                ["run_1\t0.3\t10\t8\t5\t2\t2\t0.25\t2"],
            )
            with self.assertRaises(ReportingError):
                build_analysis(directory)

    def test_summarizes_continuous_interference_degrees(self) -> None:
        rows = [
            {"run_basename": "run_1", "scan_id": "1", "classification": "low_interference", "interference_fraction": "0.05"},
            {"run_basename": "run_1", "scan_id": "2", "classification": "likely_chimeric", "interference_fraction": "0.35"},
            {"run_basename": "run_1", "scan_id": "3", "classification": "likely_chimeric", "interference_fraction": "0.90"},
            {"run_basename": "run_1", "scan_id": "4", "classification": "indeterminate", "interference_fraction": ""},
        ]
        summary = summarize_interference_distribution(rows)
        group = [row for row in summary if row["summary_level"] == "single_group"]
        self.assertEqual(5, len(group))
        self.assertEqual("3", group[0]["scorable_ms2"])
        self.assertEqual(
            {"0_to_under_10_percent", "30_to_under_50_percent", "75_to_100_percent"},
            {row["interference_band"] for row in group if row["scan_count"] == "1"},
        )
