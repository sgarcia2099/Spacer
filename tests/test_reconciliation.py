from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from spacer.bundles import discover_bundles
from spacer.reconciliation import reconcile_bundle


MZML = """<?xml version=\"1.0\" encoding=\"UTF-8\"?>
<mzML xmlns=\"http://psi.hupo.org/ms/mzml\"><run><spectrumList count=\"2\">
<spectrum id=\"controllerType=0 controllerNumber=1 scan=10\"><cvParam accession=\"MS:1000511\" name=\"ms level\" value=\"1\" /></spectrum>
<spectrum id=\"controllerType=0 controllerNumber=1 scan=11\">
  <cvParam accession=\"MS:1000511\" name=\"ms level\" value=\"2\" />
  <scanList><scan><cvParam accession=\"MS:1000016\" name=\"scan start time\" value=\"2.0\" unitName=\"minute\" /></scan></scanList>
  <precursorList><precursor><isolationWindow><cvParam accession=\"MS:1000827\" name=\"isolation window target m/z\" value=\"500.2\" /><cvParam accession=\"MS:1000828\" name=\"isolation window lower offset\" value=\"0.7\" /><cvParam accession=\"MS:1000829\" name=\"isolation window upper offset\" value=\"0.7\" /></isolationWindow><selectedIonList><selectedIon><cvParam accession=\"MS:1000744\" name=\"selected ion m/z\" value=\"500.2\" /><cvParam accession=\"MS:1000041\" name=\"charge state\" value=\"2\" /></selectedIon></selectedIonList></precursor></precursorList>
</spectrum></spectrumList></run></mzML>"""


class ReconciliationTests(unittest.TestCase):
    def test_exact_scan_metadata_reconciliation(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            stem = "run_1"
            (directory / f"{stem}.mzML").write_text(MZML, encoding="utf-8")
            (directory / f"{stem}_MSMSSpectrumInfo.txt").write_text(
                '"File ID"\t"RT in min"\t"First Scan"\t"MS Order"\t"Number of PSMs"\t"Isolation Interference in Percent"\t"Precursor mz in Da"\t"Precursor Charge"\n'
                '"F1"\t"2.0"\t"11"\t"MS2"\t"1"\t"25.0"\t"500.2"\t"2"\n',
                encoding="utf-8",
            )
            (directory / f"{stem}_PSMs.txt").write_text(
                '"First Scan"\t"File ID"\t"Annotated Sequence"\t"q-Value"\t"PEP"\n'
                '"11"\t"F1"\t"[K].PEPTIDE.[R]"\t"0.001"\t"0.002"\n',
                encoding="utf-8",
            )
            (directory / f"{stem}_PeptideGroups.txt").write_text(
                '"Annotated Sequence"\t"Number of PSMs"\n"[K].PEPTIDE.[R]"\t"1"\n',
                encoding="utf-8",
            )
            bundle = discover_bundles(input_path=None, input_dir=directory)[0]
            rows, summary = reconcile_bundle(
                bundle,
                q_value_cutoff=0.01,
                rt_tolerance_seconds=1.0,
                precursor_mz_tolerance_da=0.01,
            )
            self.assertEqual(1, len(rows))
            self.assertEqual("matched_metadata_exact_precursor", rows[0]["reconciliation_status"])
            self.assertEqual("true", rows[0]["pd_identified"])
            self.assertEqual("25.0", rows[0]["pd_isolation_interference_percent"])
            self.assertEqual("1", summary["matched_metadata_exact_precursor"])

    def test_accepts_pd_monoisotope_inside_isolation_window(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            stem = "run_1"
            (directory / f"{stem}.mzML").write_text(MZML, encoding="utf-8")
            (directory / f"{stem}_MSMSSpectrumInfo.txt").write_text(
                '"File ID"\t"RT in min"\t"First Scan"\t"MS Order"\t"Number of PSMs"\t"Precursor mz in Da"\t"Precursor Charge"\n'
                '"F1"\t"2.0"\t"11"\t"MS2"\t"0"\t"500.7"\t"2"\n',
                encoding="utf-8",
            )
            (directory / f"{stem}_PSMs.txt").write_text(
                '"First Scan"\t"File ID"\t"Annotated Sequence"\t"q-Value"\t"PEP"\n',
                encoding="utf-8",
            )
            (directory / f"{stem}_PeptideGroups.txt").write_text(
                '"Annotated Sequence"\t"Number of PSMs"\n', encoding="utf-8"
            )
            bundle = discover_bundles(input_path=None, input_dir=directory)[0]
            rows, summary = reconcile_bundle(
                bundle,
                q_value_cutoff=0.01,
                rt_tolerance_seconds=1.0,
                precursor_mz_tolerance_da=0.01,
            )
            self.assertEqual(
                "matched_metadata_precursor_within_isolation_window",
                rows[0]["reconciliation_status"],
            )
            self.assertEqual("true", rows[0]["pd_precursor_within_mzml_isolation_window"])
            self.assertEqual("1", summary["matched_metadata_precursor_within_isolation_window"])
