from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from spacer.agreement import agreement_for_bundles
from spacer.bundles import discover_bundles


class AgreementTests(unittest.TestCase):
    def test_describes_pd_context_without_modifying_score_rows(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            stem = "run_1"
            (root / f"{stem}.mzML").write_text("<mzML/>", encoding="utf-8")
            (root / f"{stem}_MSMSSpectrumInfo.txt").write_text(
                '"First Scan"\t"File ID"\t"RT in min"\t"MS Order"\t"Number of PSMs"\t"Precursor mz in Da"\t"Precursor Charge"\t"Isolation Interference in Percent"\n'
                '"2"\t"file"\t"1.0"\t"MS2"\t"1"\t"500"\t"2"\t"40"\n'
                '"3"\t"file"\t"1.1"\t"MS2"\t"0"\t"600"\t"2"\t""\n',
                encoding="utf-8",
            )
            (root / f"{stem}_PSMs.txt").write_text(
                '"First Scan"\t"File ID"\t"Annotated Sequence"\t"q-Value"\t"PEP"\n'
                '"2"\t"file"\t"PEPTIDE"\t"0.001"\t"0.01"\n',
                encoding="utf-8",
            )
            (root / f"{stem}_PeptideGroups.txt").write_text('"Annotated Sequence"\t"Number of PSMs"\n', encoding="utf-8")
            score = root / "score.tsv"
            score.write_text(
                "run_basename\tscan_id\tclassification\tinterference_fraction\ttarget_precursor_fraction\tcompetitor_count\n"
                "run_1\t2\tlikely_chimeric\t0.5\t0.5\t1\n"
                "run_1\t3\tindeterminate\t\t\t0\n",
                encoding="utf-8",
            )
            bundles = discover_bundles(input_path=None, input_dir=root)
            rows, summaries, discordant = agreement_for_bundles(
                bundles, score, q_value_cutoff=0.01, top_discordant=5
            )
            self.assertEqual("matched_finite", rows[0]["agreement_status"])
            self.assertEqual("50", rows[0]["spacer_interference_percent"])
            self.assertEqual("40", rows[0]["pd_isolation_interference_percent"])
            self.assertEqual("spacer_indeterminate", rows[1]["agreement_status"])
            self.assertEqual(4, len(summaries))
            mapped = next(row for row in summaries if row["subset"] == "pd_spectrum_mapped")
            self.assertEqual("1", mapped["scorable_ms2"])
            self.assertEqual("1", mapped["likely_chimeric"])
            self.assertEqual("1", mapped["likely_chimeric_fraction"])
            self.assertEqual("10", discordant[0]["absolute_difference_percent_points"])
