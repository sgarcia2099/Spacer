from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from spacer.conversion import ConversionError, plan_conversions


class ConversionPlanningTests(unittest.TestCase):
    def test_marks_raw_without_mzml_as_pending(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            (directory / "run_1.raw").touch()
            plan = plan_conversions(directory)[0]
            self.assertEqual("pending_conversion", plan.status)
            self.assertEqual(directory / "run_1.mzML", plan.mzml_path)

    def test_skips_matching_existing_mzml(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            (directory / "run_1.raw").touch()
            (directory / "run_1.mzml").touch()
            plan = plan_conversions(directory)[0]
            self.assertEqual("skipped_existing_mzml", plan.status)
            self.assertEqual(directory / "run_1.mzml", plan.mzml_path)

    def test_rejects_duplicate_raw_stems(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            (directory / "run_1.raw").touch()
            (directory / "RUN_1.RAW").touch()
            with self.assertRaisesRegex(ConversionError, "Ambiguous .raw"):
                plan_conversions(directory)
