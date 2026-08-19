"""Descriptive exact-scan agreement with Proteome Discoverer interference."""

from __future__ import annotations

import csv
import math
import statistics
from pathlib import Path

from .bundles import RunBundle
from .reconciliation import ReconciliationError, _read_pd_psms, _read_pd_spectra


class AgreementError(ValueError):
    """Raised when optional PD agreement inputs are invalid."""


SCORE_FIELDS = {
    "run_basename", "scan_id", "classification", "interference_fraction",
    "target_precursor_fraction", "competitor_count",
}


def _read_scores(path: Path) -> list[dict[str, str]]:
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise AgreementError(f"{path} has no header.")
            missing = SCORE_FIELDS.difference(reader.fieldnames)
            if missing:
                raise AgreementError(f"{path} is missing column(s): {', '.join(sorted(missing))}.")
            return list(reader)
    except OSError as exc:
        raise AgreementError(f"Cannot read scoring table {path}: {exc}") from exc


def _float(value: str) -> float | None:
    if not value:
        return None
    parsed = float(value)
    if not math.isfinite(parsed):
        return None
    return parsed


def _rank(values: list[float]) -> list[float]:
    ordering = sorted(enumerate(values), key=lambda item: item[1])
    ranks = [0.0] * len(values)
    start = 0
    while start < len(ordering):
        end = start + 1
        while end < len(ordering) and ordering[end][1] == ordering[start][1]:
            end += 1
        average_rank = (start + 1 + end) / 2
        for index, _ in ordering[start:end]:
            ranks[index] = average_rank
        start = end
    return ranks


def _pearson(x: list[float], y: list[float]) -> float | None:
    if len(x) < 2:
        return None
    x_mean = statistics.fmean(x)
    y_mean = statistics.fmean(y)
    numerator = sum((a - x_mean) * (b - y_mean) for a, b in zip(x, y, strict=True))
    x_ss = sum((a - x_mean) ** 2 for a in x)
    y_ss = sum((b - y_mean) ** 2 for b in y)
    if x_ss == 0 or y_ss == 0:
        return None
    return numerator / math.sqrt(x_ss * y_ss)


def _format(value: float | None) -> str:
    return "" if value is None else f"{value:.12g}"


def _summary(run: str, subset: str, rows: list[dict[str, str]], all_rows: list[dict[str, str]]) -> dict[str, str]:
    pairs = [
        (float(row["spacer_interference_percent"]), float(row["pd_isolation_interference_percent"]))
        for row in rows
        if row["agreement_status"] == "matched_finite"
    ]
    spacer = [pair[0] for pair in pairs]
    pd = [pair[1] for pair in pairs]
    differences = [a - b for a, b in pairs]
    return {
        "run_basename": run,
        "subset": subset,
        "score_rows": str(len(all_rows)),
        "matched_finite": str(len(pairs)),
        "pd_interference_missing": str(sum(row["agreement_status"] == "pd_interference_missing" for row in all_rows)),
        "spacer_indeterminate": str(sum(row["agreement_status"] == "spacer_indeterminate" for row in all_rows)),
        "scoring_without_pd_spectrum": str(sum(row["agreement_status"] == "scoring_without_pd_spectrum" for row in all_rows)),
        "pearson_r": _format(_pearson(spacer, pd)),
        "spearman_rho": _format(_pearson(_rank(spacer), _rank(pd))),
        "median_signed_difference_pp": _format(statistics.median(differences) if differences else None),
        "median_absolute_difference_pp": _format(statistics.median([abs(value) for value in differences]) if differences else None),
    }


def agreement_for_bundles(
    bundles: list[RunBundle], scoring_path: Path, *, q_value_cutoff: float, top_discordant: int
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    """Return exact-scan PD context and descriptive agreement only."""
    if not 0 <= q_value_cutoff <= 1 or top_discordant < 1:
        raise AgreementError("q-value cutoff must be 0-1 and top discordant count must be positive.")
    score_rows = _read_scores(scoring_path)
    bundle_map = {bundle.run_basename: bundle for bundle in bundles}
    grouped: dict[str, dict[int, dict[str, str]]] = {}
    for row in score_rows:
        try:
            scan = int(row["scan_id"])
        except ValueError as exc:
            raise AgreementError(f"Invalid scoring scan ID: {row['scan_id']!r}") from exc
        if row["run_basename"] not in bundle_map:
            raise AgreementError(f"Scoring table references unknown run {row['run_basename']!r}.")
        if scan in grouped.setdefault(row["run_basename"], {}):
            raise AgreementError(f"Duplicate scoring scan {scan} in run {row['run_basename']!r}.")
        grouped[row["run_basename"]][scan] = row
    if set(grouped) != set(bundle_map):
        raise AgreementError("Scoring table and input bundles have different run basenames.")
    agreement_rows: list[dict[str, str]] = []
    summaries: list[dict[str, str]] = []
    discordant: list[dict[str, str]] = []
    for run, bundle in bundle_map.items():
        try:
            pd_spectra = _read_pd_spectra(bundle.msms_spectrum_info)
            psm = _read_pd_psms(bundle.psms, q_value_cutoff)
        except ReconciliationError as exc:
            raise AgreementError(str(exc)) from exc
        run_rows: list[dict[str, str]] = []
        for scan, score in grouped[run].items():
            pd = pd_spectra.get(scan)
            spacer_fraction = _float(score["interference_fraction"])
            pd_value = pd.isolation_interference_percent if pd is not None else None
            if pd is None:
                status = "scoring_without_pd_spectrum"
            elif spacer_fraction is None:
                status = "spacer_indeterminate"
            elif pd_value is None:
                status = "pd_interference_missing"
            else:
                status = "matched_finite"
            difference = 100 * spacer_fraction - pd_value if status == "matched_finite" and spacer_fraction is not None and pd_value is not None else None
            psm_summary = psm.get(scan)
            row = {
                "run_basename": run,
                "scan_id": str(scan),
                "agreement_status": status,
                "score_classification": score["classification"],
                "spacer_interference_percent": _format(100 * spacer_fraction if spacer_fraction is not None else None),
                "pd_isolation_interference_percent": _format(pd_value),
                "signed_difference_percent_points": _format(difference),
                "absolute_difference_percent_points": _format(abs(difference) if difference is not None else None),
                "pd_identified": str(psm_summary.identified if psm_summary else False).lower(),
                "pd_best_q_value": _format(psm_summary.best_q_value if psm_summary else None),
                "pd_raw_psm_count": str(psm_summary.raw_psm_count if psm_summary else 0),
                "competitor_count": score["competitor_count"],
            }
            run_rows.append(row)
        agreement_rows.extend(run_rows)
        summaries.append(_summary(run, "all_matched", run_rows, run_rows))
        summaries.append(_summary(run, "identified", [row for row in run_rows if row["pd_identified"] == "true"], run_rows))
        summaries.append(_summary(run, "unidentified", [row for row in run_rows if row["pd_identified"] == "false"], run_rows))
        matched = [row for row in run_rows if row["agreement_status"] == "matched_finite"]
        discordant.extend(sorted(matched, key=lambda row: float(row["absolute_difference_percent_points"]), reverse=True)[:top_discordant])
    return agreement_rows, summaries, sorted(
        discordant, key=lambda row: float(row["absolute_difference_percent_points"]), reverse=True
    )

