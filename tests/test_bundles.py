from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from spacer.bundles import BundleValidationError, discover_bundles


MSMS_HEADER = '"First Scan"\t"MS Order"\t"Isolation Interference in Percent"\n'
PSM_HEADER = '"First Scan"\t"File ID"\t"Annotated Sequence"\t"q-Value"\t"PEP"\n'
PEPTIDE_HEADER = '"Annotated Sequence"\t"Number of PSMs"\n'


def _write_bundle(directory: Path, stem: str, *, missing: str | None = None) -> None:
    (directory / f"{stem}.raw").touch()
    files = {
        "msms_spectrum_info": (f"{stem}_MSMSSpectrumInfo.txt", MSMS_HEADER),
        "psms": (f"{stem}_PSMs.txt", PSM_HEADER),
        "peptide_groups": (f"{stem}_PeptideGroups.txt", PEPTIDE_HEADER),
    }
    for kind, (name, header) in files.items():
        if kind != missing:
            (directory / name).write_text(header, encoding="utf-8")


class BundleDiscoveryTests(unittest.TestCase):
    def test_discovers_one_complete_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            _write_bundle(directory, "run_1")
            bundles = discover_bundles(input_path=None, input_dir=directory)
            self.assertEqual(1, len(bundles))
            self.assertEqual("run_1", bundles[0].run_basename)
            self.assertTrue(bundles[0].has_pd_isolation_interference)

    def test_prefers_existing_mzml_over_raw(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            _write_bundle(directory, "run_1")
            (directory / "run_1.mzML").touch()
            bundle = discover_bundles(input_path=None, input_dir=directory)[0]
            self.assertEqual("mzml", bundle.data_format)

    def test_rejects_incomplete_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            _write_bundle(directory, "run_1", missing="psms")
            with self.assertRaisesRegex(BundleValidationError, "missing result export"):
                discover_bundles(input_path=None, input_dir=directory)

    def test_rejects_missing_required_header(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            _write_bundle(directory, "run_1")
            (directory / "run_1_PSMs.txt").write_text('"First Scan"\n', encoding="utf-8")
            with self.assertRaisesRegex(BundleValidationError, "missing required column"):
                discover_bundles(input_path=None, input_dir=directory)
