"""Concise descriptive Analysis reports for one DDA technical-replicate group."""

from __future__ import annotations

import csv
import statistics
from pathlib import Path


class ReportingError(ValueError):
    """Raised when completed scoring outputs cannot form a descriptive report."""


RUN_FIELDS = {
    "run_basename",
    "total_ms2",
    "scorable_ms2",
    "low_interference",
    "likely_chimeric",
    "indeterminate",
    "likely_chimeric_fraction",
    "interference_threshold",
}
SENSITIVITY_FIELDS = RUN_FIELDS | {"minimum_competitor_prior_ms1_detections"}
INTERFERENCE_BANDS = (
    ("0_to_under_10_percent", 0.00, 0.10),
    ("10_to_under_30_percent", 0.10, 0.30),
    ("30_to_under_50_percent", 0.30, 0.50),
    ("50_to_under_75_percent", 0.50, 0.75),
    ("75_to_100_percent", 0.75, 1.00),
)


def _read_tsv(path: Path, required: set[str]) -> list[dict[str, str]]:
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise ReportingError(f"{path} has no tab-separated header.")
            missing = required.difference(reader.fieldnames)
            if missing:
                raise ReportingError(f"{path} is missing column(s): {', '.join(sorted(missing))}.")
            rows = list(reader)
    except OSError as exc:
        raise ReportingError(f"Cannot read {path}: {exc}") from exc
    if not rows:
        raise ReportingError(f"{path} has no data rows.")
    return rows


def integer(row: dict[str, str], field: str) -> int:
    try:
        return int(row[field])
    except (KeyError, ValueError) as exc:
        raise ReportingError(f"Invalid integer {field!r} in run {row.get('run_basename', '?')!r}.") from exc


def _fraction(row: dict[str, str], field: str) -> float:
    try:
        value = float(row[field])
    except (KeyError, ValueError) as exc:
        raise ReportingError(f"Invalid fraction {field!r} in run {row.get('run_basename', '?')!r}.") from exc
    if not 0 <= value <= 1:
        raise ReportingError(f"Fraction {field!r} is outside 0-1 in run {row.get('run_basename', '?')!r}.")
    return value


def _validate_run_row(row: dict[str, str]) -> None:
    total = integer(row, "total_ms2")
    scorable = integer(row, "scorable_ms2")
    low = integer(row, "low_interference")
    likely = integer(row, "likely_chimeric")
    indeterminate = integer(row, "indeterminate")
    if total != scorable + indeterminate or scorable != low + likely:
        raise ReportingError(f"Scoring arithmetic is inconsistent for run {row['run_basename']!r}.")
    expected = likely / scorable if scorable else None
    if expected is not None and abs(_fraction(row, "likely_chimeric_fraction") - expected) > 1e-12:
        raise ReportingError(f"Likely-chimeric fraction is inconsistent for run {row['run_basename']!r}.")


def _format(value: float | None) -> str:
    return "" if value is None else f"{value:.12g}"


def _percentile(values: list[float], fraction: float) -> float:
    """Return an interpolated percentile without a numerical dependency."""
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (position - lower)


def summarize_interference_distribution(
    score_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    """Summarize the continuous independent interference estimate by run and group."""
    grouped: dict[str, list[float]] = {}
    for row in score_rows:
        if row.get("classification") == "indeterminate":
            continue
        try:
            value = float(row["interference_fraction"])
        except (KeyError, ValueError) as exc:
            raise ReportingError(
                f"Scorable scan {row.get('scan_id', '?')!r} has no valid interference fraction."
            ) from exc
        if not 0 <= value <= 1:
            raise ReportingError(
                f"Interference fraction is outside 0-1 for scan {row.get('scan_id', '?')!r}."
            )
        grouped.setdefault(row["run_basename"], []).append(value)
    if not grouped:
        raise ReportingError("No scorable interference fractions are available.")
    grouped["single_group"] = [
        value for run, values in grouped.items() if run != "single_group" for value in values
    ]
    rows: list[dict[str, str]] = []
    for run, values in grouped.items():
        for label, lower, upper in INTERFERENCE_BANDS:
            count = sum(
                lower <= value < upper or (upper == 1 and value == 1)
                for value in values
            )
            rows.append(
                {
                    "summary_level": "single_group" if run == "single_group" else "run",
                    "run_basename": run,
                    "scorable_ms2": str(len(values)),
                    "interference_band": label,
                    "lower_inclusive_percent": str(100 * lower),
                    "upper_inclusive_percent": str(100 * upper),
                    "scan_count": str(count),
                    "fraction_of_scorable_ms2": _format(count / len(values)),
                    "median_interference_fraction": _format(_percentile(values, 0.5)),
                    "p10_interference_fraction": _format(_percentile(values, 0.1)),
                    "p25_interference_fraction": _format(_percentile(values, 0.25)),
                    "p75_interference_fraction": _format(_percentile(values, 0.75)),
                    "p90_interference_fraction": _format(_percentile(values, 0.9)),
                }
            )
    return rows


def build_analysis(
    scoring_dir: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]], str]:
    """Build per-run and single-group descriptive summaries from score outputs."""
    run_rows = _read_tsv(scoring_dir / "run_summary.tsv", RUN_FIELDS)
    sensitivity_rows = _read_tsv(
        scoring_dir / "sensitivity_summary.tsv", SENSITIVITY_FIELDS
    )
    run_names = [row["run_basename"] for row in run_rows]
    if len(set(run_names)) != len(run_names):
        raise ReportingError("run_summary.tsv contains duplicate run basenames.")
    for row in run_rows:
        _validate_run_row(row)

    analysis_rows: list[dict[str, str]] = []
    for row in run_rows:
        total = integer(row, "total_ms2")
        indeterminate = integer(row, "indeterminate")
        analysis_rows.append(
            {
                "summary_level": "run",
                "run_basename": row["run_basename"],
                "run_count": "1",
                "total_ms2": str(total),
                "scorable_ms2": row["scorable_ms2"],
                "low_interference": row["low_interference"],
                "likely_chimeric": row["likely_chimeric"],
                "indeterminate": row["indeterminate"],
                "likely_chimeric_fraction": row["likely_chimeric_fraction"],
                "indeterminate_fraction": _format(indeterminate / total if total else None),
                "replicate_likely_fraction_median": "",
                "replicate_likely_fraction_min": "",
                "replicate_likely_fraction_max": "",
                "interference_threshold": row["interference_threshold"],
            }
        )
    total = sum(integer(row, "total_ms2") for row in run_rows)
    scorable = sum(integer(row, "scorable_ms2") for row in run_rows)
    low = sum(integer(row, "low_interference") for row in run_rows)
    likely = sum(integer(row, "likely_chimeric") for row in run_rows)
    indeterminate = sum(integer(row, "indeterminate") for row in run_rows)
    replicate_fractions = [_fraction(row, "likely_chimeric_fraction") for row in run_rows]
    thresholds = {row["interference_threshold"] for row in run_rows}
    if len(thresholds) != 1:
        raise ReportingError("run_summary.tsv has inconsistent interference thresholds.")
    analysis_rows.append(
        {
            "summary_level": "single_group",
            "run_basename": "single_group",
            "run_count": str(len(run_rows)),
            "total_ms2": str(total),
            "scorable_ms2": str(scorable),
            "low_interference": str(low),
            "likely_chimeric": str(likely),
            "indeterminate": str(indeterminate),
            "likely_chimeric_fraction": _format(likely / scorable if scorable else None),
            "indeterminate_fraction": _format(indeterminate / total if total else None),
            "replicate_likely_fraction_median": _format(statistics.median(replicate_fractions)),
            "replicate_likely_fraction_min": _format(min(replicate_fractions)),
            "replicate_likely_fraction_max": _format(max(replicate_fractions)),
            "interference_threshold": thresholds.pop(),
        }
    )

    sensitivity_by_threshold: dict[str, list[dict[str, str]]] = {}
    for row in sensitivity_rows:
        _validate_run_row(row)
        sensitivity_by_threshold.setdefault(row["interference_threshold"], []).append(row)
    sensitivity_summary: list[dict[str, str]] = []
    for threshold, rows in sorted(sensitivity_by_threshold.items(), key=lambda item: float(item[0])):
        if {row["run_basename"] for row in rows} != set(run_names):
            raise ReportingError(f"Sensitivity threshold {threshold} does not cover every run.")
        fractions = [_fraction(row, "likely_chimeric_fraction") for row in rows]
        summary_total = sum(integer(row, "total_ms2") for row in rows)
        summary_scorable = sum(integer(row, "scorable_ms2") for row in rows)
        summary_likely = sum(integer(row, "likely_chimeric") for row in rows)
        summary_indeterminate = sum(integer(row, "indeterminate") for row in rows)
        sensitivity_summary.append(
            {
                "interference_threshold": threshold,
                "run_count": str(len(rows)),
                "total_ms2": str(summary_total),
                "scorable_ms2": str(summary_scorable),
                "likely_chimeric": str(summary_likely),
                "indeterminate": str(summary_indeterminate),
                "pooled_likely_chimeric_fraction": _format(
                    summary_likely / summary_scorable if summary_scorable else None
                ),
                "replicate_likely_fraction_median": _format(statistics.median(fractions)),
                "replicate_likely_fraction_min": _format(min(fractions)),
                "replicate_likely_fraction_max": _format(max(fractions)),
            }
        )
    group = analysis_rows[-1]
    report = "\n".join(
        [
            "# Spacer Analysis report",
            "",
            "Single-group descriptive summary; no between-method comparison or hypothesis test was run.",
            "",
            f"- Runs: {group['run_count']}",
            f"- Total MS2: {group['total_ms2']}",
            f"- Scorable MS2: {group['scorable_ms2']}",
            f"- Likely chimeric: {group['likely_chimeric_fraction']} (pooled)",
            f"- Replicate likely-chimeric range: {group['replicate_likely_fraction_min']} to {group['replicate_likely_fraction_max']}",
            f"- Indeterminate: {group['indeterminate_fraction']} (pooled)",
            f"- Primary interference threshold: {group['interference_threshold']}",
            "",
        ]
    )
    return analysis_rows, sensitivity_summary, report


def write_tsv(
    path: Path, rows: list[dict[str, str]], *, fieldnames: list[str] | None = None
) -> None:
    """Write TSV rows, retaining a useful header when a result is empty."""
    if not rows and fieldnames is None:
        raise ReportingError(f"Cannot write {path}: there are no rows or field names.")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames or list(rows[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)
