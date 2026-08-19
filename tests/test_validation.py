from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from spacer.bundles import discover_bundles
from spacer.validation import validate_completed_analysis


class ValidationModeTests(unittest.TestCase):
    def test_validates_completed_single_run_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            stem = "run_1"
            (root / f"{stem}.mzML").write_text(
                '<mzML><run><spectrumList><spectrum id="scan=1"><cvParam accession="MS:1000511" value="1" /></spectrum><spectrum id="scan=2"><cvParam accession="MS:1000511" value="2" /></spectrum></spectrumList></run></mzML>',
                encoding="utf-8",
            )
            (root / f"{stem}_MSMSSpectrumInfo.txt").write_text('"First Scan"\t"MS Order"\n"2"\t"MS2"\n', encoding="utf-8")
            (root / f"{stem}_PSMs.txt").write_text('"First Scan"\t"File ID"\t"Annotated Sequence"\t"q-Value"\t"PEP"\n', encoding="utf-8")
            (root / f"{stem}_PeptideGroups.txt").write_text('"Annotated Sequence"\t"Number of PSMs"\n', encoding="utf-8")
            scoring = root / "score"
            scoring.mkdir()
            (scoring / "run_summary.tsv").write_text(
                "run_basename\ttotal_ms2\tscorable_ms2\tlow_interference\tlikely_chimeric\tindeterminate\tlikely_chimeric_fraction\tinterference_threshold\n"
                "run_1\t1\t1\t1\t0\t0\t0.0\t0.3\n",
                encoding="utf-8",
            )
            (scoring / "sensitivity_summary.tsv").write_text(
                "run_basename\tinterference_threshold\ttotal_ms2\tscorable_ms2\tlow_interference\tlikely_chimeric\tindeterminate\tlikely_chimeric_fraction\tminimum_competitor_prior_ms1_detections\n"
                "run_1\t0.3\t1\t1\t1\t0\t0\t0.0\t2\n",
                encoding="utf-8",
            )
            (scoring / "ms2_chimericity.tsv").write_text(
                "run_basename\tscan_id\tprecursor_mz\treported_charge\tisolation_lower_mz\tisolation_upper_mz\ttarget_precursor_fraction\tinterference_fraction\tstrongest_competitor_fraction\tcompetitor_count\ttarget_prior_ms1_detections\tcompetitor_prior_ms1_detections\tclassification\tindeterminate_reason\n"
                "run_1\t2\t500\t2\t499\t501\t1.0\t0.0\t0.0\t0\t1\t0\tlow_interference\t\n",
                encoding="utf-8",
            )
            bundles = discover_bundles(input_path=None, input_dir=root)
            checks, report = validate_completed_analysis(bundles, scoring, None)
            self.assertEqual(5, len(checks))
            self.assertIn("Status: PASS", report)
