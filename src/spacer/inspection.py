"""Opt-in candidate selection and coordinate export for manual inspection."""

from __future__ import annotations

import csv
import heapq
import math
import xml.etree.ElementTree as ET
from pathlib import Path

from .bundles import RunBundle
from .reconciliation import _local_name
from .scoring import (
    _envelope_indices,
    _ms_level,
    _peak_arrays,
    _scan_number,
    _spectrum_params,
)


class InspectionError(ValueError):
    """Raised when inspection inputs are incomplete or incompatible."""


SCORE_FIELDS = {
    "run_basename",
    "scan_id",
    "precursor_mz",
    "reported_charge",
    "isolation_lower_mz",
    "isolation_upper_mz",
    "interference_fraction",
    "classification",
    "indeterminate_reason",
}


def _number(value: str, name: str) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        raise InspectionError(f"Score row has an invalid {name}: {value!r}.") from exc
    if not math.isfinite(parsed):
        raise InspectionError(f"Score row has a non-finite {name}: {value!r}.")
    return parsed


def read_score_rows(path: Path) -> list[dict[str, str]]:
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise InspectionError(f"Scoring table has no header: {path}")
            missing = SCORE_FIELDS.difference(reader.fieldnames)
            if missing:
                raise InspectionError(
                    f"Scoring table is missing required column(s): {', '.join(sorted(missing))}."
                )
            rows = list(reader)
    except OSError as exc:
        raise InspectionError(f"Cannot read scoring table {path}: {exc}") from exc
    if not rows:
        raise InspectionError(f"Scoring table has no MS2 rows: {path}")
    return rows


def read_discordant_rows(path: Path) -> list[dict[str, str]]:
    required = {
        "run_basename", "scan_id", "agreement_status",
        "absolute_difference_percent_points",
    }
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise InspectionError(f"Agreement table has no header: {path}")
            missing = required.difference(reader.fieldnames)
            if missing:
                raise InspectionError(
                    f"Agreement table is missing required column(s): {', '.join(sorted(missing))}."
                )
            return list(reader)
    except OSError as exc:
        raise InspectionError(f"Cannot read agreement table {path}: {exc}") from exc


def _score_value(row: dict[str, str]) -> float:
    value = row["interference_fraction"]
    return _number(value, "interference_fraction") if value else float("nan")


def select_candidates(
    rows: list[dict[str, str]], *, interference_threshold: float, per_category: int,
    balance_runs: bool = False, discordant_rows: list[dict[str, str]] | None = None,
) -> list[dict[str, str]]:
    """Select auditable representatives without rendering a plot."""
    if per_category < 1:
        raise InspectionError("Candidate count per category must be positive.")
    with_score = [row for row in rows if row["interference_fraction"]]
    groups: tuple[tuple[str, str, list[dict[str, str]], object, bool], ...] = (
        (
            "high_interference",
            "highest independent interference among likely-chimeric scans",
            [row for row in with_score if row["classification"] == "likely_chimeric"],
            _score_value,
            True,
        ),
        (
            "low_interference",
            "lowest independent interference among low-interference scans",
            [row for row in with_score if row["classification"] == "low_interference"],
            _score_value,
            False,
        ),
        (
            "threshold_adjacent",
            f"nearest independent interference to configured threshold {interference_threshold:g}",
            with_score,
            lambda row: abs(_score_value(row) - interference_threshold),
            False,
        ),
        (
            "indeterminate",
            "indeterminate score requiring manual review",
            [row for row in rows if row["classification"] == "indeterminate"],
            lambda row: (row["indeterminate_reason"], row["run_basename"], int(row["scan_id"])),
            False,
        ),
    )
    selected: list[dict[str, str]] = []
    for category, rationale, candidates, key, largest in groups:
        candidates = _select_ranked(
            candidates, limit=per_category, key=key, largest=largest, balance_runs=balance_runs
        )
        for rank, row in enumerate(candidates, start=1):
            selected.append(
                {
                    "candidate_category": category,
                    "category_rank": str(rank),
                    "candidate_scope": "balanced_across_runs" if balance_runs else "global",
                    "selection_rationale": rationale,
                    "pd_agreement_status": "",
                    "pd_absolute_difference_percent_points": "",
                    **row,
                }
            )
    if discordant_rows:
        score_by_key = {(row["run_basename"], row["scan_id"]): row for row in rows}
        pool = [
                row for row in discordant_rows
                if row["agreement_status"] == "matched_finite"
                and row["absolute_difference_percent_points"]
        ]
        ranked = _select_ranked(
            pool,
            limit=per_category,
            key=lambda row: _number(
                row["absolute_difference_percent_points"],
                "PD absolute difference percent points",
            ),
            largest=True,
            balance_runs=balance_runs,
        )
        for rank, agreement in enumerate(ranked, start=1):
            score = score_by_key.get((agreement["run_basename"], agreement["scan_id"]))
            if score is None:
                raise InspectionError(
                    f"Agreement candidate is absent from scoring table: {agreement['run_basename']} scan {agreement['scan_id']}."
                )
            selected.append(
                {
                    "candidate_category": "pd_discordant",
                    "category_rank": str(rank),
                    "candidate_scope": "balanced_across_runs" if balance_runs else "global",
                    "selection_rationale": "largest exact-scan Spacer versus PD interference difference",
                    "pd_agreement_status": agreement["agreement_status"],
                    "pd_absolute_difference_percent_points": agreement["absolute_difference_percent_points"],
                    **score,
                }
            )
    return selected


def _select_ranked(
    candidates: list[dict[str, str]], *, limit: int, key: object, largest: bool, balance_runs: bool
) -> list[dict[str, str]]:
    if not balance_runs:
        selector = heapq.nlargest if largest else heapq.nsmallest
        return selector(limit, candidates, key=key)
    """Ensure each run contributes its best candidate before global repeats."""
    selected: list[dict[str, str]] = []
    selected_ids: set[int] = set()
    selector = heapq.nlargest if largest else heapq.nsmallest
    for run in sorted({row["run_basename"] for row in candidates}):
        run_rows = [row for row in candidates if row["run_basename"] == run]
        if run_rows and len(selected) < limit:
            chosen = selector(1, run_rows, key=key)[0]
            selected.append(chosen)
            selected_ids.add(id(chosen))
    remaining = [row for row in candidates if id(row) not in selected_ids]
    selected.extend(selector(limit - len(selected), remaining, key=key))
    return selected


def _array_rows(
    *, run_basename: str, scan_id: int, preceding_ms1_scan_id: int | None,
    mz: tuple[float, ...], intensity: tuple[float, ...], lower: float, upper: float,
    margin: float, target_indices: set[int],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for index, (peak_mz, peak_intensity) in enumerate(zip(mz, intensity, strict=True)):
        if lower - margin <= peak_mz <= upper + margin:
            rows.append(
                {
                    "run_basename": run_basename,
                    "ms2_scan_id": str(scan_id),
                    "preceding_ms1_scan_id": str(preceding_ms1_scan_id or ""),
                    "mz": str(peak_mz),
                    "intensity": str(peak_intensity),
                    "in_isolation_window": "true" if lower <= peak_mz <= upper else "false",
                    "is_target_envelope": "true" if index in target_indices else "false",
                    "is_transmitted_target_signal": (
                        "true" if index in target_indices and lower <= peak_mz <= upper else "false"
                    ),
                }
            )
    return rows


def export_scan_coordinates(
    bundle: RunBundle, *, score_row: dict[str, str], output_dir: Path,
    ppm_tolerance: float, max_isotopes: int, window_margin_mz: float,
) -> dict[str, Path]:
    """Export source coordinates for one explicitly selected MS2 scan."""
    if bundle.data_format != "mzml":
        raise InspectionError(f"Run {bundle.run_basename} requires mzML; run spacer convert first.")
    try:
        requested_scan = int(score_row["scan_id"])
        lower = _number(score_row["isolation_lower_mz"], "isolation_lower_mz")
        upper = _number(score_row["isolation_upper_mz"], "isolation_upper_mz")
        target_mz = _number(score_row["precursor_mz"], "precursor_mz")
        charge = int(score_row["reported_charge"])
    except (ValueError, InspectionError) as exc:
        raise InspectionError(f"Scan {score_row['scan_id']} cannot be inspected: {exc}") from exc
    prior_ms1: tuple[int | None, tuple[float, ...], tuple[float, ...]] | None = None
    ms2_arrays: tuple[tuple[float, ...], tuple[float, ...]] | None = None
    found = False
    try:
        for _, spectrum in ET.iterparse(bundle.data_path, events=("end",)):
            if _local_name(spectrum.tag) != "spectrum":
                continue
            params = _spectrum_params(spectrum)
            level = _ms_level(params)
            scan_id = _scan_number(spectrum.attrib.get("id", ""))
            if level == 1:
                mz, intensity = _peak_arrays(spectrum)
                prior_ms1 = (scan_id, mz, intensity)
            elif level == 2 and scan_id == requested_scan:
                ms2_arrays = _peak_arrays(spectrum)
                found = True
                spectrum.clear()
                break
            spectrum.clear()
    except (ET.ParseError, OSError, ValueError) as exc:
        raise InspectionError(f"Cannot read {bundle.data_path} for scan {requested_scan}: {exc}") from exc
    if not found or ms2_arrays is None:
        raise InspectionError(f"MS2 scan {requested_scan} was not found in {bundle.data_path}.")
    if prior_ms1 is None:
        raise InspectionError(f"MS2 scan {requested_scan} has no preceding MS1 survey scan.")
    preceding_ms1_scan_id, ms1_mz, ms1_intensity = prior_ms1
    target_index = min(range(len(ms1_mz)), key=lambda index: abs(ms1_mz[index] - target_mz), default=None)
    target_indices: set[int] = set()
    if target_index is not None:
        target_indices = set(_envelope_indices(ms1_mz, target_index, charge, ppm_tolerance, max_isotopes))
    ms1_rows = _array_rows(
        run_basename=bundle.run_basename,
        scan_id=requested_scan,
        preceding_ms1_scan_id=preceding_ms1_scan_id,
        mz=ms1_mz,
        intensity=ms1_intensity,
        lower=lower,
        upper=upper,
        margin=window_margin_mz,
        target_indices=target_indices,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{bundle.run_basename}_scan_{requested_scan}"
    ms1_path = output_dir / f"{stem}_ms1_window_coordinates.tsv"
    ms2_path = output_dir / f"{stem}_ms2_coordinates.tsv"
    metadata_path = output_dir / f"{stem}_inspection_metadata.tsv"
    _write_tsv(
        ms1_path,
        ms1_rows,
        fieldnames=[
            "run_basename", "ms2_scan_id", "preceding_ms1_scan_id", "mz", "intensity",
            "in_isolation_window", "is_target_envelope",
            "is_transmitted_target_signal",
        ],
    )
    _write_tsv(
        ms2_path,
        [
            {
                "run_basename": bundle.run_basename,
                "ms2_scan_id": str(requested_scan),
                "mz": str(mz),
                "intensity": str(intensity),
            }
            for mz, intensity in zip(*ms2_arrays, strict=True)
        ],
        fieldnames=["run_basename", "ms2_scan_id", "mz", "intensity"],
    )
    metadata = {
        "run_basename": bundle.run_basename,
        "ms2_scan_id": str(requested_scan),
        "preceding_ms1_scan_id": str(preceding_ms1_scan_id or ""),
        "isolation_lower_mz": str(lower),
        "isolation_upper_mz": str(upper),
        "window_margin_mz": str(window_margin_mz),
        "score_classification": score_row["classification"],
        "score_interference_fraction": score_row["interference_fraction"],
        "score_indeterminate_reason": score_row["indeterminate_reason"],
        "plot_generated": "false",
    }
    _write_tsv(metadata_path, [{"field": key, "value": value} for key, value in metadata.items()])
    return {"ms1": ms1_path, "ms2": ms2_path, "metadata": metadata_path}


def collect_scan_plot_data(
    bundle: RunBundle, *, score_rows: list[dict[str, str]], ppm_tolerance: float,
    max_isotopes: int, window_margin_mz: float,
) -> list[dict[str, object]]:
    """Collect coordinate arrays for selected scans in one pass through one mzML file."""
    if bundle.data_format != "mzml":
        raise InspectionError(f"Run {bundle.run_basename} requires mzML; run spacer convert first.")
    requested: dict[int, dict[str, str]] = {}
    for row in score_rows:
        try:
            scan = int(row["scan_id"])
        except (KeyError, ValueError) as exc:
            raise InspectionError("Selected report plot has an invalid scan ID.") from exc
        if scan in requested:
            raise InspectionError(f"Selected report plot duplicates scan {scan} in {bundle.run_basename}.")
        requested[scan] = row
    prior_ms1: tuple[int | None, tuple[float, ...], tuple[float, ...]] | None = None
    plots: list[dict[str, object]] = []
    try:
        for _, spectrum in ET.iterparse(bundle.data_path, events=("end",)):
            if _local_name(spectrum.tag) != "spectrum":
                continue
            params = _spectrum_params(spectrum)
            level = _ms_level(params)
            scan_id = _scan_number(spectrum.attrib.get("id", ""))
            if level == 1:
                mz, intensity = _peak_arrays(spectrum)
                prior_ms1 = (scan_id, mz, intensity)
            elif level == 2 and scan_id in requested:
                if prior_ms1 is None:
                    raise InspectionError(f"MS2 scan {scan_id} has no preceding MS1 survey scan.")
                row = requested.pop(scan_id)
                lower = _number(row["isolation_lower_mz"], "isolation_lower_mz")
                upper = _number(row["isolation_upper_mz"], "isolation_upper_mz")
                target_mz = _number(row["precursor_mz"], "precursor_mz")
                try:
                    charge = int(row["reported_charge"])
                except ValueError as exc:
                    raise InspectionError(f"Scan {scan_id} has an invalid reported charge.") from exc
                preceding_ms1_scan_id, ms1_mz, ms1_intensity = prior_ms1
                target_index = min(
                    range(len(ms1_mz)),
                    key=lambda index: abs(ms1_mz[index] - target_mz),
                    default=None,
                )
                target_indices: set[int] = set()
                if target_index is not None:
                    target_indices = set(
                        _envelope_indices(
                            ms1_mz, target_index, charge, ppm_tolerance, max_isotopes
                        )
                    )
                ms1_rows = _array_rows(
                    run_basename=bundle.run_basename,
                    scan_id=scan_id,
                    preceding_ms1_scan_id=preceding_ms1_scan_id,
                    mz=ms1_mz,
                    intensity=ms1_intensity,
                    lower=lower,
                    upper=upper,
                    margin=window_margin_mz,
                    target_indices=target_indices,
                )
                ms2_mz, ms2_intensity = _peak_arrays(spectrum)
                plots.append(
                    {
                        "run_basename": bundle.run_basename,
                        "scan_id": str(scan_id),
                        "preceding_ms1_scan_id": str(preceding_ms1_scan_id or ""),
                        "isolation_lower_mz": str(lower),
                        "isolation_upper_mz": str(upper),
                        "classification": row["classification"],
                        "interference_fraction": row["interference_fraction"],
                        "ms1": ms1_rows,
                        "ms2": [
                            {"mz": str(mz), "intensity": str(intensity)}
                            for mz, intensity in zip(ms2_mz, ms2_intensity, strict=True)
                        ],
                    }
                )
            spectrum.clear()
    except (ET.ParseError, OSError, ValueError) as exc:
        raise InspectionError(f"Cannot collect report plots from {bundle.data_path}: {exc}") from exc
    if requested:
        missing = ", ".join(str(scan) for scan in sorted(requested))
        raise InspectionError(f"Selected MS2 scan(s) not found in {bundle.data_path}: {missing}.")
    return sorted(plots, key=lambda row: int(str(row["scan_id"])))


def render_scan_plot(paths: dict[str, Path]) -> dict[str, Path]:
    """Render an explicitly requested two-panel MS1/MS2 inspection plot."""
    try:
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
    except ImportError as exc:
        raise InspectionError(
            "Plot rendering requires matplotlib. Install the optional plot extra first."
        ) from exc
    ms1_rows = _read_tsv(paths["ms1"])
    ms2_rows = _read_tsv(paths["ms2"])
    metadata = {row["field"]: row["value"] for row in _read_tsv(paths["metadata"])}
    if not ms1_rows or not ms2_rows:
        raise InspectionError("Cannot plot an inspection export with no MS1 or MS2 coordinate rows.")
    lower = _number(metadata["isolation_lower_mz"], "isolation_lower_mz")
    upper = _number(metadata["isolation_upper_mz"], "isolation_upper_mz")
    scan = metadata["ms2_scan_id"]
    interference = 100 * _number(metadata["score_interference_fraction"], "score_interference_fraction")
    figure, (ms1_axis, ms2_axis) = plt.subplots(1, 2, figsize=(14, 5.5), constrained_layout=True)
    ms1_axis.axvspan(lower, upper, color="#f4a261", alpha=0.18)
    for row in ms1_rows:
        mz = _number(row["mz"], "MS1 m/z")
        intensity = _number(row["intensity"], "MS1 intensity")
        in_window = row["in_isolation_window"] == "true"
        target = row["is_target_envelope"] == "true"
        if target and in_window:
            color, style = "#00a86b", "-"
        elif target:
            color, style = "#2a9d8f", "--"
        elif in_window:
            color, style = "#303030", "-"
        else:
            color, style = "#b5b5b5", "-"
        ms1_axis.vlines(mz, 1, max(intensity, 1), color=color, linestyle=style, linewidth=1.3)
    ms1_axis.set_yscale("log")
    ms1_axis.set_xlabel("m/z")
    ms1_axis.set_ylabel("MS1 intensity (log scale)")
    ms1_axis.set_title(f"MS1 isolation window | scan {scan}\ninterference {interference:.1f}%")
    ms1_axis.legend(
        handles=[
            Line2D([0], [0], color="#303030", label="Other peak in window"),
            Line2D([0], [0], color="#b5b5b5", label="Peak outside window"),
            Line2D([0], [0], color="#00a86b", label="Target envelope transmitted"),
            Line2D([0], [0], color="#2a9d8f", linestyle="--", label="Target envelope outside nominal window"),
        ],
        fontsize=8,
        loc="upper left",
    )
    ms2_mz = [_number(row["mz"], "MS2 m/z") for row in ms2_rows]
    ms2_intensity = [_number(row["intensity"], "MS2 intensity") for row in ms2_rows]
    maximum = max(ms2_intensity)
    relative = [100 * value / maximum if maximum else 0 for value in ms2_intensity]
    ms2_axis.vlines(ms2_mz, 0, relative, color="#305f9f", linewidth=0.8)
    ms2_axis.set_ylim(0, 105)
    ms2_axis.set_xlabel("m/z")
    ms2_axis.set_ylabel("MS2 relative intensity (%)")
    ms2_axis.set_title(f"MS2 spectrum | {metadata['score_classification']}")
    stem = paths["metadata"].name.removesuffix("_inspection_metadata.tsv")
    png_path = paths["metadata"].with_name(f"{stem}_inspection_plot.png")
    svg_path = paths["metadata"].with_name(f"{stem}_inspection_plot.svg")
    figure.savefig(png_path, dpi=200)
    figure.savefig(svg_path)
    plt.close(figure)
    return {"png": png_path, "svg": svg_path}


def _read_tsv(path: Path) -> list[dict[str, str]]:
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))
    except OSError as exc:
        raise InspectionError(f"Cannot read inspection output {path}: {exc}") from exc


def _write_tsv(
    path: Path, rows: list[dict[str, str]], *, fieldnames: list[str] | None = None
) -> None:
    if not rows and fieldnames is None:
        raise InspectionError(f"No rows were available for {path.name}.")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames or list(rows[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)
