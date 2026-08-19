"""Detailed structural and result-consistency checks for Validation mode."""

from __future__ import annotations

import csv
import math
from pathlib import Path

from .bundles import RunBundle
from .conversion import ConversionError, inspect_mzml
from .inspection import read_score_rows
from .reporting import ReportingError, build_analysis, integer, write_tsv


class ValidationError(ValueError):
    """Raised when a Validation-mode check fails."""


def _read_tsv(path: Path, required: set[str]) -> list[dict[str, str]]:
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise ValidationError(f"{path} has no header.")
            missing = required.difference(reader.fieldnames)
            if missing:
                raise ValidationError(f"{path} is missing column(s): {', '.join(sorted(missing))}.")
            return list(reader)
    except OSError as exc:
        raise ValidationError(f"Cannot read {path}: {exc}") from exc


def _check(run: str, component: str, detail: str) -> dict[str, str]:
    return {
        "run_basename": run,
        "component": component,
        "severity": "info",
        "status": "pass",
        "detail": detail,
    }


def validate_completed_analysis(
    bundles: list[RunBundle], scoring_dir: Path, reconciliation_dir: Path | None,
    agreement_dir: Path | None = None,
) -> tuple[list[dict[str, str]], str]:
    """Validate existing scoring outputs against bundle and mzML structure."""
    try:
        analysis_rows, sensitivity_rows, _ = build_analysis(scoring_dir)
    except ReportingError as exc:
        raise ValidationError(str(exc)) from exc
    run_rows = [row for row in analysis_rows if row["summary_level"] == "run"]
    scoring_summary = {row["run_basename"]: row for row in run_rows}
    bundle_map = {bundle.run_basename: bundle for bundle in bundles}
    if set(scoring_summary) != set(bundle_map):
        raise ValidationError("Scoring summaries and discovered input bundles have different run basenames.")
    checks: list[dict[str, str]] = []
    score_rows = read_score_rows(scoring_dir / "ms2_chimericity.tsv")
    by_run: dict[str, list[dict[str, str]]] = {}
    for row in score_rows:
        by_run.setdefault(row["run_basename"], []).append(row)
    if set(by_run) != set(bundle_map):
        raise ValidationError("MS2 scoring table and discovered input bundles have different run basenames.")
    for run, bundle in bundle_map.items():
        if bundle.data_format != "mzml":
            raise ValidationError(f"Validation requires mzML for run {run}; run spacer convert first.")
        try:
            mzml = inspect_mzml(bundle.data_path)
        except ConversionError as exc:
            raise ValidationError(str(exc)) from exc
        rows = by_run[run]
        scan_ids = [row["scan_id"] for row in rows]
        if any(not scan for scan in scan_ids) or len(scan_ids) != len(set(scan_ids)):
            raise ValidationError(f"Scoring table has missing or duplicate scan IDs for run {run}.")
        summary = scoring_summary[run]
        if len(rows) != mzml.ms2_count or len(rows) != integer(summary, "total_ms2"):
            raise ValidationError(
                f"MS2 count disagreement for {run}: mzML={mzml.ms2_count}, scoring rows={len(rows)}, summary={summary['total_ms2']}."
            )
        classified = {"low_interference": 0, "likely_chimeric": 0, "indeterminate": 0}
        for row in rows:
            classification = row["classification"]
            if classification not in classified:
                raise ValidationError(f"Unknown scoring classification {classification!r} in run {run}.")
            classified[classification] += 1
            for field in ("target_precursor_fraction", "interference_fraction", "strongest_competitor_fraction"):
                if row[field]:
                    value = float(row[field])
                    if not math.isfinite(value) or not 0 <= value <= 1:
                        raise ValidationError(f"Invalid {field} in run {run}, scan {row['scan_id']}.")
            if classification == "likely_chimeric" and row["indeterminate_reason"]:
                raise ValidationError(f"Likely-chimeric row has an indeterminate reason in run {run}.")
            if classification == "indeterminate" and not row["indeterminate_reason"]:
                raise ValidationError(f"Indeterminate row lacks a reason in run {run}.")
        for name, value in classified.items():
            if value != integer(summary, name):
                raise ValidationError(f"Summary count mismatch for {run}, {name}.")
        checks.extend(
            [
                _check(run, "bundle", "Complete Single-mode bundle discovered."),
                _check(run, "mzml_structure", f"MS1={mzml.ms1_count}; MS2={mzml.ms2_count}; total={mzml.spectrum_count}."),
                _check(run, "score_rows", f"{len(rows)} unique MS2 score rows agree with mzML and run summary."),
                _check(run, "score_values", "Classifications, fractions, and indeterminate reasons are internally consistent."),
            ]
        )

    thresholds = {row["interference_threshold"] for row in sensitivity_rows}
    checks.append(
        _check(
            "single_group",
            "sensitivity_grid",
            f"Validated {len(thresholds)} prespecified threshold(s): {', '.join(sorted(thresholds, key=float))}.",
        )
    )
    if reconciliation_dir is not None:
        reconciliation_rows = _read_tsv(
            reconciliation_dir / "scan_reconciliation_summary.tsv",
            {
                "run_basename", "mzml_ms2_count", "pd_ms2_count",
                "matched_metadata_exact_precursor",
                "matched_metadata_precursor_within_isolation_window",
                "matched_scan_metadata_mismatch_or_missing",
                "mzml_ms2_without_pd_spectrum", "pd_spectrum_without_mzml_ms2",
            },
        )
        reconciliation = {row["run_basename"]: row for row in reconciliation_rows}
        if set(reconciliation) != set(bundle_map):
            raise ValidationError("Reconciliation summaries and discovered input bundles have different run basenames.")
        for run, row in reconciliation.items():
            union_count = sum(
                int(row[field])
                for field in (
                    "matched_metadata_exact_precursor",
                    "matched_metadata_precursor_within_isolation_window",
                    "matched_scan_metadata_mismatch_or_missing",
                    "mzml_ms2_without_pd_spectrum",
                    "pd_spectrum_without_mzml_ms2",
                )
            )
            if union_count < int(row["mzml_ms2_count"]):
                raise ValidationError(f"Reconciliation union count is smaller than mzML MS2 count for {run}.")
            checks.append(
                _check(
                    run,
                    "reconciliation",
                    f"mzML MS2={row['mzml_ms2_count']}; PD MS2={row['pd_ms2_count']}; union scan keys={union_count}.",
                )
            )
    if agreement_dir is not None:
        agreement_rows = _read_tsv(
            agreement_dir / "pd_agreement_summary.tsv",
            {
                "run_basename", "subset", "score_rows", "matched_finite",
                "pd_interference_missing", "spacer_indeterminate",
                "scoring_without_pd_spectrum", "pearson_r", "spearman_rho",
            },
        )
        all_matched = {
            row["run_basename"]: row for row in agreement_rows if row["subset"] == "all_matched"
        }
        if set(all_matched) != set(bundle_map):
            raise ValidationError("Agreement outputs need exactly one all_matched summary for every run.")
        for run, row in all_matched.items():
            if int(row["score_rows"]) != integer(scoring_summary[run], "total_ms2"):
                raise ValidationError(f"PD agreement score-row count mismatch for {run}.")
            finite = int(row["matched_finite"])
            if finite < 0 or finite > int(row["score_rows"]):
                raise ValidationError(f"PD agreement matched count is invalid for {run}.")
            for field in ("pearson_r", "spearman_rho"):
                if row[field] and not -1 <= float(row[field]) <= 1:
                    raise ValidationError(f"PD agreement {field} is outside -1 to 1 for {run}.")
            checks.append(
                _check(
                    run,
                    "pd_agreement",
                    f"Matched finite={row['matched_finite']}; PD missing={row['pd_interference_missing']}; Spacer indeterminate={row['spacer_indeterminate']}.",
                )
            )
    report = "\n".join(
        [
            "# Spacer Validation report",
            "",
            f"Status: PASS ({len(checks)} checks)",
            "",
            "This report validates declared structure and arithmetic. It does not establish biological truth or recalibrate Spacer scores from Proteome Discoverer.",
            "",
        ]
    )
    return checks, report


def write_validation(output_dir: Path, checks: list[dict[str, str]], report: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(output_dir / "validation_checks.tsv", checks)
    (output_dir / "validation_report.md").write_text(report, encoding="utf-8")
