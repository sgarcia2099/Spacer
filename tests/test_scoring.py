from __future__ import annotations

import base64
import struct
import tempfile
import unittest
import zlib
from pathlib import Path

from spacer.bundles import discover_bundles
from spacer.scoring import score_bundle, sensitivity_summary


def _array(kind_accession: str, values: list[float]) -> str:
    payload = zlib.compress(struct.pack(f"<{len(values)}d", *values))
    return (
        f'<binaryDataArray><cvParam accession="{kind_accession}" />'
        '<cvParam accession="MS:1000523" /><cvParam accession="MS:1000574" />'
        f"<binary>{base64.b64encode(payload).decode()}</binary></binaryDataArray>"
    )


class ScoringTests(unittest.TestCase):
    def test_scores_competing_isotope_envelope_without_pd_input(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            stem = "run_1"
            # Target: z=2 at 500.0000/500.5017. A stronger, distinct z=2
            # envelope falls inside the same isolation window.
            mz = [499.2, 499.7017, 500.0, 500.2034, 500.5017, 500.7051]
            intensity = [200.0, 200.0, 100.0, 200.0, 50.0, 200.0]
            mzml = f'''<mzML><run><spectrumList>
<spectrum id="scan=1"><cvParam accession="MS:1000511" value="1" /><binaryDataArrayList>{_array("MS:1000514", mz)}{_array("MS:1000515", intensity)}</binaryDataArrayList></spectrum>
<spectrum id="scan=2"><cvParam accession="MS:1000511" value="2" /><precursorList><precursor><isolationWindow><cvParam accession="MS:1000827" value="500.0" /><cvParam accession="MS:1000828" value="1.0" /><cvParam accession="MS:1000829" value="1.0" /></isolationWindow><selectedIonList><selectedIon><cvParam accession="MS:1000041" value="2" /></selectedIon></selectedIonList></precursor></precursorList></spectrum>
</spectrumList></run></mzML>'''
            (directory / f"{stem}.mzML").write_text(mzml, encoding="utf-8")
            (directory / f"{stem}_MSMSSpectrumInfo.txt").write_text('"First Scan"\t"MS Order"\n"2"\t"MS2"\n', encoding="utf-8")
            (directory / f"{stem}_PSMs.txt").write_text('"First Scan"\t"File ID"\t"Annotated Sequence"\t"q-Value"\t"PEP"\n', encoding="utf-8")
            (directory / f"{stem}_PeptideGroups.txt").write_text('"Annotated Sequence"\t"Number of PSMs"\n', encoding="utf-8")
            bundle = discover_bundles(input_path=None, input_dir=directory)[0]
            rows, summary = score_bundle(
                bundle, interference_threshold=0.30, ppm_tolerance=20, max_isotopes=3,
                previous_ms1_scans=1, minimum_competitor_prior_ms1_detections=1,
            )
            self.assertEqual(1, len(rows))
            self.assertEqual("likely_chimeric", rows[0]["classification"])
            self.assertGreater(float(rows[0]["interference_fraction"]), 0.30)
            self.assertGreaterEqual(int(rows[0]["competitor_count"]), 1)
            self.assertEqual("1", summary["total_ms2"])
            self.assertEqual("1", rows[0]["competitor_prior_ms1_detections"])

            strict_rows, _ = score_bundle(
                bundle, interference_threshold=0.30, ppm_tolerance=20, max_isotopes=3,
                previous_ms1_scans=1, minimum_competitor_prior_ms1_detections=2,
            )
            self.assertEqual("indeterminate", strict_rows[0]["classification"])
            self.assertEqual(
                "competitor_insufficient_coelution_support",
                strict_rows[0]["indeterminate_reason"],
            )
            sensitivity = sensitivity_summary(
                rows, thresholds=(0.20, 0.30, 0.40),
                minimum_competitor_prior_ms1_detections=1,
            )
            self.assertEqual(["0.2", "0.3", "0.4"], [row["interference_threshold"] for row in sensitivity])

    def test_marks_ms2_before_any_ms1_as_indeterminate(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            stem = "run_2"
            mzml = '''<mzML><run><spectrumList>
<spectrum id="scan=1"><cvParam accession="MS:1000511" value="2" /><precursorList><precursor><isolationWindow><cvParam accession="MS:1000827" value="500.0" /><cvParam accession="MS:1000828" value="1.0" /><cvParam accession="MS:1000829" value="1.0" /></isolationWindow><selectedIonList><selectedIon><cvParam accession="MS:1000041" value="2" /></selectedIon></selectedIonList></precursor></precursorList></spectrum>
</spectrumList></run></mzML>'''
            (directory / f"{stem}.mzML").write_text(mzml, encoding="utf-8")
            (directory / f"{stem}_MSMSSpectrumInfo.txt").write_text('"First Scan"\t"MS Order"\n"1"\t"MS2"\n', encoding="utf-8")
            (directory / f"{stem}_PSMs.txt").write_text('"First Scan"\t"File ID"\t"Annotated Sequence"\t"q-Value"\t"PEP"\n', encoding="utf-8")
            (directory / f"{stem}_PeptideGroups.txt").write_text('"Annotated Sequence"\t"Number of PSMs"\n', encoding="utf-8")
            bundle = discover_bundles(input_path=None, input_dir=directory)[0]
            rows, summary = score_bundle(
                bundle, interference_threshold=0.30, ppm_tolerance=20, max_isotopes=3,
                previous_ms1_scans=1, minimum_competitor_prior_ms1_detections=1,
            )
            self.assertEqual("indeterminate", rows[0]["classification"])
            self.assertEqual("no_preceding_ms1", rows[0]["indeterminate_reason"])
            self.assertEqual("0", summary["scorable_ms2"])
