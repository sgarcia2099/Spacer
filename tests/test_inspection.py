from __future__ import annotations

import base64
import importlib.util
import struct
import tempfile
import unittest
import zlib
from pathlib import Path

from spacer.bundles import discover_bundles
from spacer.inspection import (
    collect_scan_plot_data,
    export_scan_coordinates,
    render_scan_plot,
    select_candidates,
)


def _array(kind_accession: str, values: list[float]) -> str:
    payload = zlib.compress(struct.pack(f"<{len(values)}d", *values))
    return (
        f'<binaryDataArray><cvParam accession="{kind_accession}" />'
        '<cvParam accession="MS:1000523" /><cvParam accession="MS:1000574" />'
        f"<binary>{base64.b64encode(payload).decode()}</binary></binaryDataArray>"
    )


class InspectionTests(unittest.TestCase):
    def test_selects_each_required_candidate_category(self) -> None:
        def row(scan_id: int, classification: str, interference: str, reason: str = "") -> dict[str, str]:
            return {
                "run_basename": "run_1", "scan_id": str(scan_id), "precursor_mz": "500",
                "reported_charge": "2", "isolation_lower_mz": "499",
                "isolation_upper_mz": "501", "interference_fraction": interference,
                "classification": classification, "indeterminate_reason": reason,
            }

        candidates = select_candidates(
            [
                row(1, "likely_chimeric", "0.9"),
                row(2, "low_interference", "0.01"),
                row(3, "low_interference", "0.31"),
                row(4, "indeterminate", "", "missing_or_invalid_precursor_charge"),
            ],
            interference_threshold=0.30,
            per_category=1,
        )
        self.assertEqual(
            {"high_interference", "low_interference", "threshold_adjacent", "indeterminate"},
            {candidate["candidate_category"] for candidate in candidates},
        )

    def test_balanced_candidates_represent_each_run_before_repeats(self) -> None:
        rows = [
            {
                "run_basename": run, "scan_id": str(index), "precursor_mz": "500",
                "reported_charge": "2", "isolation_lower_mz": "499",
                "isolation_upper_mz": "501", "interference_fraction": str(0.9 - index / 100),
                "classification": "likely_chimeric", "indeterminate_reason": "",
            }
            for index, run in enumerate(("run_1", "run_2", "run_3"), start=1)
        ]
        candidates = select_candidates(
            rows, interference_threshold=0.30, per_category=3, balance_runs=True
        )
        high = [row for row in candidates if row["candidate_category"] == "high_interference"]
        self.assertEqual({"run_1", "run_2", "run_3"}, {row["run_basename"] for row in high})
        self.assertTrue(all(row["candidate_scope"] == "balanced_across_runs" for row in high))

    def test_adds_pd_discordant_candidate_only_from_exact_finite_agreement(self) -> None:
        score = {
            "run_basename": "run_1", "scan_id": "2", "precursor_mz": "500",
            "reported_charge": "2", "isolation_lower_mz": "499",
            "isolation_upper_mz": "501", "interference_fraction": "0.9",
            "classification": "likely_chimeric", "indeterminate_reason": "",
        }
        candidates = select_candidates(
            [score],
            interference_threshold=0.30,
            per_category=1,
            discordant_rows=[
                {
                    "run_basename": "run_1", "scan_id": "2",
                    "agreement_status": "matched_finite",
                    "absolute_difference_percent_points": "80",
                },
                {
                    "run_basename": "run_1", "scan_id": "3",
                    "agreement_status": "spacer_indeterminate",
                    "absolute_difference_percent_points": "",
                },
            ],
        )
        discordant = [row for row in candidates if row["candidate_category"] == "pd_discordant"]
        self.assertEqual(1, len(discordant))
        self.assertEqual("80", discordant[0]["pd_absolute_difference_percent_points"])

    def test_exports_coordinates_without_plot(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            stem = "run_1"
            mzml = f'''<mzML><run><spectrumList>
<spectrum id="scan=1"><cvParam accession="MS:1000511" value="1" /><binaryDataArrayList>{_array("MS:1000514", [498.8, 500.0, 500.5017, 501.2])}{_array("MS:1000515", [5.0, 100.0, 50.0, 3.0])}</binaryDataArrayList></spectrum>
<spectrum id="scan=2"><cvParam accession="MS:1000511" value="2" /><binaryDataArrayList>{_array("MS:1000514", [100.0, 200.0])}{_array("MS:1000515", [10.0, 20.0])}</binaryDataArrayList></spectrum>
</spectrumList></run></mzML>'''
            (directory / f"{stem}.mzML").write_text(mzml, encoding="utf-8")
            (directory / f"{stem}_MSMSSpectrumInfo.txt").write_text('"First Scan"\t"MS Order"\n"2"\t"MS2"\n', encoding="utf-8")
            (directory / f"{stem}_PSMs.txt").write_text('"First Scan"\t"File ID"\t"Annotated Sequence"\t"q-Value"\t"PEP"\n', encoding="utf-8")
            (directory / f"{stem}_PeptideGroups.txt").write_text('"Annotated Sequence"\t"Number of PSMs"\n', encoding="utf-8")
            bundle = discover_bundles(input_path=None, input_dir=directory)[0]
            paths = export_scan_coordinates(
                bundle,
                score_row={
                    "scan_id": "2", "precursor_mz": "500", "reported_charge": "2",
                    "isolation_lower_mz": "499", "isolation_upper_mz": "501",
                    "classification": "likely_chimeric", "interference_fraction": "0.5",
                    "indeterminate_reason": "",
                },
                output_dir=directory / "inspection",
                ppm_tolerance=20,
                max_isotopes=3,
                window_margin_mz=0.5,
            )
            self.assertTrue(paths["ms1"].is_file())
            self.assertTrue(paths["ms2"].is_file())
            self.assertIn("plot_generated\tfalse", paths["metadata"].read_text(encoding="utf-8"))
            self.assertIn("is_target_envelope", paths["ms1"].read_text(encoding="utf-8"))
            self.assertIn("is_transmitted_target_signal", paths["ms1"].read_text(encoding="utf-8"))
            collected = collect_scan_plot_data(
                bundle,
                score_rows=[
                    {
                        "scan_id": "2", "precursor_mz": "500", "reported_charge": "2",
                        "isolation_lower_mz": "499", "isolation_upper_mz": "501",
                        "classification": "likely_chimeric", "interference_fraction": "0.5",
                    }
                ],
                ppm_tolerance=20,
                max_isotopes=3,
                window_margin_mz=0.5,
            )
            self.assertEqual("2", collected[0]["scan_id"])
            self.assertTrue(collected[0]["ms1"])
            self.assertTrue(collected[0]["ms2"])
            if importlib.util.find_spec("matplotlib") is not None:
                plots = render_scan_plot(paths)
                self.assertTrue(plots["png"].is_file())
                self.assertTrue(plots["svg"].is_file())
