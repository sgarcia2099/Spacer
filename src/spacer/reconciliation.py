"""Exact scan reconciliation between mzML and Proteome Discoverer exports."""

from __future__ import annotations

import csv
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

from .bundles import RunBundle


class ReconciliationError(ValueError):
    """Raised when required scan metadata cannot be parsed or reconciled."""


SCAN_PATTERN = re.compile(r"(?:^|\s)scan=(\d+)(?:\s|$)")


@dataclass(frozen=True)
class MzmlMs2:
    scan_id: int
    native_id: str
    rt_seconds: float | None
    precursor_mz: float | None
    precursor_charge: int | None
    isolation_target_mz: float | None
    isolation_lower_offset: float | None
    isolation_upper_offset: float | None


@dataclass(frozen=True)
class PdSpectrum:
    scan_id: int
    file_id: str
    rt_minutes: float | None
    precursor_mz: float | None
    precursor_charge: int | None
    psm_count: int | None
    isolation_interference_percent: float | None


@dataclass(frozen=True)
class PdPsmSummary:
    raw_psm_count: int
    best_q_value: float | None
    identified: bool


def _local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _param(params: dict[str, ET.Element], accession: str, name: str) -> ET.Element | None:
    """Get a cvParam without relying on Element truthiness."""

    found = params.get(accession)
    return found if found is not None else params.get(name)


def _as_float(value: str | None, context: str) -> float | None:
    if value is None or not value.strip():
        return None
    try:
        return float(value)
    except ValueError as exc:
        raise ReconciliationError(f"Invalid numeric value for {context}: {value!r}") from exc


def _as_int(value: str | None, context: str) -> int | None:
    parsed = _as_float(value, context)
    if parsed is None:
        return None
    if not parsed.is_integer():
        raise ReconciliationError(f"Expected integer value for {context}: {value!r}")
    return int(parsed)


def _header_reader(path: Path) -> tuple[object, csv.DictReader]:
    try:
        handle = path.open("r", encoding="utf-8-sig", newline="")
    except OSError as exc:
        raise ReconciliationError(f"Cannot read {path}: {exc}") from exc
    reader = csv.DictReader(handle, delimiter="\t")
    if reader.fieldnames is None:
        handle.close()
        raise ReconciliationError(f"{path} has no tab-separated header.")
    reader.fieldnames = [name.strip() for name in reader.fieldnames]
    return handle, reader


def _read_mzml_ms2(path: Path) -> dict[int, MzmlMs2]:
    scans: dict[int, MzmlMs2] = {}
    try:
        for _, spectrum in ET.iterparse(path, events=("end",)):
            if _local_name(spectrum.tag) != "spectrum":
                continue
            params = {
                child.attrib.get("accession", child.attrib.get("name", "")): child
                for child in spectrum.iter()
                if _local_name(child.tag) == "cvParam"
            }
            level = _param(params, "MS:1000511", "ms level")
            if level is None or level.attrib.get("value") != "2":
                spectrum.clear()
                continue
            native_id = spectrum.attrib.get("id", "")
            match = SCAN_PATTERN.search(native_id)
            if match is None:
                raise ReconciliationError(
                    f"Cannot extract scan=<number> from mzML MS2 native ID: {native_id!r}"
                )
            scan_id = int(match.group(1))
            if scan_id in scans:
                raise ReconciliationError(f"Duplicate MS2 scan number {scan_id} in {path}")

            rt = _param(params, "MS:1000016", "scan start time")
            rt_value = _as_float(rt.attrib.get("value") if rt is not None else None, "mzML scan start time")
            if rt_value is not None and rt is not None:
                unit = (rt.attrib.get("unitName") or rt.attrib.get("unitAccession") or "").casefold()
                if unit in {"minute", "min", "uo:0000031"}:
                    rt_value *= 60
                elif unit not in {"second", "sec", "s", "uo:0000010", ""}:
                    raise ReconciliationError(f"Unsupported mzML retention-time unit: {unit!r}")

            def param_value(accession: str, name: str, context: str) -> float | None:
                param = _param(params, accession, name)
                return _as_float(param.attrib.get("value") if param is not None else None, context)

            precursor_charge_value = param_value("MS:1000041", "charge state", "mzML precursor charge")
            precursor_charge = (
                int(precursor_charge_value) if precursor_charge_value is not None else None
            )
            scans[scan_id] = MzmlMs2(
                scan_id=scan_id,
                native_id=native_id,
                rt_seconds=rt_value,
                precursor_mz=param_value("MS:1000744", "selected ion m/z", "mzML precursor m/z"),
                precursor_charge=precursor_charge,
                isolation_target_mz=param_value("MS:1000827", "isolation window target m/z", "mzML isolation target m/z"),
                isolation_lower_offset=param_value("MS:1000828", "isolation window lower offset", "mzML isolation lower offset"),
                isolation_upper_offset=param_value("MS:1000829", "isolation window upper offset", "mzML isolation upper offset"),
            )
            spectrum.clear()
    except (ET.ParseError, OSError) as exc:
        raise ReconciliationError(f"Cannot parse mzML {path}: {exc}") from exc
    if not scans:
        raise ReconciliationError(f"No MS2 scans found in {path}")
    return scans


def _read_pd_spectra(path: Path) -> dict[int, PdSpectrum]:
    required = {
        "First Scan",
        "File ID",
        "RT in min",
        "MS Order",
        "Number of PSMs",
        "Precursor mz in Da",
        "Precursor Charge",
    }
    handle, reader = _header_reader(path)
    missing = sorted(required - set(reader.fieldnames or []))
    if missing:
        handle.close()
        raise ReconciliationError(f"{path} missing required column(s): {', '.join(missing)}")
    records: dict[int, PdSpectrum] = {}
    try:
        for row_number, row in enumerate(reader, start=2):
            if (row.get("MS Order") or "").strip().upper() != "MS2":
                continue
            scan_id = _as_int(row.get("First Scan"), f"{path} row {row_number} First Scan")
            if scan_id is None:
                raise ReconciliationError(f"{path} row {row_number} has no First Scan")
            if scan_id in records:
                raise ReconciliationError(f"Duplicate MS2 First Scan {scan_id} in {path}")
            interference = _as_float(
                row.get("Isolation Interference in Percent"),
                f"{path} row {row_number} Isolation Interference in Percent",
            )
            if interference is not None and not 0 <= interference <= 100:
                raise ReconciliationError(
                    f"{path} row {row_number} Isolation Interference in Percent outside 0–100: {interference}"
                )
            records[scan_id] = PdSpectrum(
                scan_id=scan_id,
                file_id=(row.get("File ID") or "").strip(),
                rt_minutes=_as_float(row.get("RT in min"), f"{path} row {row_number} RT in min"),
                precursor_mz=_as_float(row.get("Precursor mz in Da"), f"{path} row {row_number} Precursor mz in Da"),
                precursor_charge=_as_int(row.get("Precursor Charge"), f"{path} row {row_number} Precursor Charge"),
                psm_count=_as_int(row.get("Number of PSMs"), f"{path} row {row_number} Number of PSMs"),
                isolation_interference_percent=interference,
            )
    finally:
        handle.close()
    return records


def _read_pd_psms(path: Path, q_value_cutoff: float) -> dict[int, PdPsmSummary]:
    required = {"First Scan", "File ID", "q-Value", "PEP"}
    handle, reader = _header_reader(path)
    missing = sorted(required - set(reader.fieldnames or []))
    if missing:
        handle.close()
        raise ReconciliationError(f"{path} missing required column(s): {', '.join(missing)}")
    grouped: dict[int, list[float | None]] = {}
    try:
        for row_number, row in enumerate(reader, start=2):
            scan_id = _as_int(row.get("First Scan"), f"{path} row {row_number} First Scan")
            if scan_id is None:
                raise ReconciliationError(f"{path} row {row_number} has no First Scan")
            q_value = _as_float(row.get("q-Value"), f"{path} row {row_number} q-Value")
            grouped.setdefault(scan_id, []).append(q_value)
    finally:
        handle.close()
    output: dict[int, PdPsmSummary] = {}
    for scan_id, q_values in grouped.items():
        finite_q_values = [value for value in q_values if value is not None]
        best_q = min(finite_q_values) if finite_q_values else None
        output[scan_id] = PdPsmSummary(
            raw_psm_count=len(q_values),
            best_q_value=best_q,
            identified=best_q is not None and best_q <= q_value_cutoff,
        )
    return output


def _read_peptide_group_count(path: Path) -> int:
    required = {"Annotated Sequence", "Number of PSMs"}
    handle, reader = _header_reader(path)
    missing = sorted(required - set(reader.fieldnames or []))
    if missing:
        handle.close()
        raise ReconciliationError(f"{path} missing required column(s): {', '.join(missing)}")
    try:
        return sum(1 for _ in reader)
    finally:
        handle.close()


def _string(value: object | None) -> str:
    return "" if value is None else str(value)


def reconcile_bundle(
    bundle: RunBundle,
    *,
    q_value_cutoff: float,
    rt_tolerance_seconds: float,
    precursor_mz_tolerance_da: float,
) -> tuple[list[dict[str, str]], dict[str, str]]:
    """Create an exact scan reconciliation table for one converted run bundle."""

    if bundle.data_format != "mzml":
        raise ReconciliationError(
            f"Run {bundle.run_basename} has no selected mzML input. Run `spacer convert` first."
        )
    mzml_scans = _read_mzml_ms2(bundle.data_path)
    pd_scans = _read_pd_spectra(bundle.msms_spectrum_info)
    psm_summaries = _read_pd_psms(bundle.psms, q_value_cutoff)
    peptide_group_count = _read_peptide_group_count(bundle.peptide_groups)

    rows: list[dict[str, str]] = []
    status_counts: dict[str, int] = {}
    for scan_id in sorted(set(mzml_scans) | set(pd_scans)):
        mzml = mzml_scans.get(scan_id)
        pd = pd_scans.get(scan_id)
        psm = psm_summaries.get(scan_id)
        if mzml is None:
            status = "pd_spectrum_without_mzml_ms2"
        elif pd is None:
            status = "mzml_ms2_without_pd_spectrum"
        else:
            rt_difference = (
                abs(mzml.rt_seconds - pd.rt_minutes * 60)
                if mzml.rt_seconds is not None and pd.rt_minutes is not None
                else None
            )
            mz_difference = (
                abs(mzml.precursor_mz - pd.precursor_mz)
                if mzml.precursor_mz is not None and pd.precursor_mz is not None
                else None
            )
            charge_match = (
                mzml.precursor_charge == pd.precursor_charge
                if mzml.precursor_charge is not None and pd.precursor_charge is not None
                else None
            )
            precursor_within_isolation_window = (
                mzml.isolation_target_mz - mzml.isolation_lower_offset - precursor_mz_tolerance_da
                <= pd.precursor_mz
                <= mzml.isolation_target_mz + mzml.isolation_upper_offset + precursor_mz_tolerance_da
                if mzml.isolation_target_mz is not None
                and mzml.isolation_lower_offset is not None
                and mzml.isolation_upper_offset is not None
                and pd.precursor_mz is not None
                else None
            )
            if (
                rt_difference is not None
                and mz_difference is not None
                and charge_match is not None
                and rt_difference <= rt_tolerance_seconds
                and mz_difference <= precursor_mz_tolerance_da
                and charge_match
            ):
                status = "matched_metadata_exact_precursor"
            elif (
                rt_difference is not None
                and rt_difference <= rt_tolerance_seconds
                and charge_match
                and precursor_within_isolation_window is True
            ):
                status = "matched_metadata_precursor_within_isolation_window"
            else:
                status = "matched_scan_metadata_mismatch_or_missing"
        status_counts[status] = status_counts.get(status, 0) + 1
        rt_difference = (
            abs(mzml.rt_seconds - pd.rt_minutes * 60)
            if mzml is not None and pd is not None and mzml.rt_seconds is not None and pd.rt_minutes is not None
            else None
        )
        mz_difference = (
            abs(mzml.precursor_mz - pd.precursor_mz)
            if mzml is not None and pd is not None and mzml.precursor_mz is not None and pd.precursor_mz is not None
            else None
        )
        charge_match = (
            mzml.precursor_charge == pd.precursor_charge
            if mzml is not None and pd is not None and mzml.precursor_charge is not None and pd.precursor_charge is not None
            else None
        )
        precursor_within_isolation_window = (
            mzml.isolation_target_mz - mzml.isolation_lower_offset - precursor_mz_tolerance_da
            <= pd.precursor_mz
            <= mzml.isolation_target_mz + mzml.isolation_upper_offset + precursor_mz_tolerance_da
            if mzml is not None
            and pd is not None
            and mzml.isolation_target_mz is not None
            and mzml.isolation_lower_offset is not None
            and mzml.isolation_upper_offset is not None
            and pd.precursor_mz is not None
            else None
        )
        rows.append(
            {
                "run_basename": bundle.run_basename,
                "scan_id": str(scan_id),
                "reconciliation_status": status,
                "mzml_native_id": _string(mzml.native_id if mzml else None),
                "mzml_rt_seconds": _string(mzml.rt_seconds if mzml else None),
                "pd_rt_minutes": _string(pd.rt_minutes if pd else None),
                "rt_difference_seconds": _string(rt_difference),
                "mzml_precursor_mz": _string(mzml.precursor_mz if mzml else None),
                "pd_precursor_mz": _string(pd.precursor_mz if pd else None),
                "precursor_mz_difference_da": _string(mz_difference),
                "pd_precursor_within_mzml_isolation_window": _string(
                    precursor_within_isolation_window
                ).lower(),
                "mzml_precursor_charge": _string(mzml.precursor_charge if mzml else None),
                "pd_precursor_charge": _string(pd.precursor_charge if pd else None),
                "precursor_charge_match": _string(charge_match).lower(),
                "pd_file_id": _string(pd.file_id if pd else None),
                "pd_spectrum_psm_count": _string(pd.psm_count if pd else None),
                "pd_raw_psm_count": _string(psm.raw_psm_count if psm else 0),
                "pd_best_q_value": _string(psm.best_q_value if psm else None),
                "pd_identified": _string(psm.identified if psm else False).lower(),
                "pd_isolation_interference_percent": _string(
                    pd.isolation_interference_percent if pd else None
                ),
                "mzml_isolation_target_mz": _string(mzml.isolation_target_mz if mzml else None),
                "mzml_isolation_lower_offset": _string(mzml.isolation_lower_offset if mzml else None),
                "mzml_isolation_upper_offset": _string(mzml.isolation_upper_offset if mzml else None),
            }
        )
    summary = {
        "run_basename": bundle.run_basename,
        "mzml_ms2_count": str(len(mzml_scans)),
        "pd_ms2_count": str(len(pd_scans)),
        "pd_psm_scan_count": str(len(psm_summaries)),
        "pd_peptide_group_count": str(peptide_group_count),
        "q_value_cutoff": str(q_value_cutoff),
        "rt_tolerance_seconds": str(rt_tolerance_seconds),
        "precursor_mz_tolerance_da": str(precursor_mz_tolerance_da),
        "matched_metadata_exact_precursor": str(status_counts.get("matched_metadata_exact_precursor", 0)),
        "matched_metadata_precursor_within_isolation_window": str(
            status_counts.get("matched_metadata_precursor_within_isolation_window", 0)
        ),
        "matched_scan_metadata_mismatch_or_missing": str(status_counts.get("matched_scan_metadata_mismatch_or_missing", 0)),
        "mzml_ms2_without_pd_spectrum": str(status_counts.get("mzml_ms2_without_pd_spectrum", 0)),
        "pd_spectrum_without_mzml_ms2": str(status_counts.get("pd_spectrum_without_mzml_ms2", 0)),
    }
    return rows, summary
