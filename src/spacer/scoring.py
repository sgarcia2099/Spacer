"""Independent, conservative MS1 precursor-interference scoring for DDA MS2."""

from __future__ import annotations

import base64
import bisect
import struct
import xml.etree.ElementTree as ET
import zlib
from collections import deque
from dataclasses import dataclass

from .bundles import RunBundle
from .reconciliation import _local_name, _param


class ScoringError(ValueError):
    """Raised when mzML spectra cannot support independent precursor scoring."""


ISOTOPE_MASS = 1.0033548378


@dataclass(frozen=True)
class Ms1Scan:
    scan_id: int | None
    rt_seconds: float | None
    mz: tuple[float, ...]
    intensity: tuple[float, ...]


def _float(value: str | None) -> float | None:
    if value is None or not value.strip():
        return None
    return float(value)


def _scan_number(native_id: str) -> int | None:
    marker = "scan="
    if marker not in native_id:
        return None
    suffix = native_id.split(marker, 1)[1].split()[0]
    try:
        return int(suffix)
    except ValueError:
        return None


def _spectrum_params(spectrum: ET.Element) -> dict[str, ET.Element]:
    return {
        child.attrib.get("accession", child.attrib.get("name", "")): child
        for child in spectrum.iter()
        if _local_name(child.tag) == "cvParam"
    }


def _decode_array(binary_data_array: ET.Element) -> tuple[str | None, tuple[float, ...]]:
    accessions = {
        child.attrib.get("accession", "")
        for child in binary_data_array
        if _local_name(child.tag) == "cvParam"
    }
    if "MS:1000514" in accessions:
        kind = "mz"
    elif "MS:1000515" in accessions:
        kind = "intensity"
    else:
        return None, ()
    binary = next((child for child in binary_data_array if _local_name(child.tag) == "binary"), None)
    if binary is None or not (binary.text or "").strip():
        return kind, ()
    payload = base64.b64decode((binary.text or "").strip())
    if "MS:1000574" in accessions:
        payload = zlib.decompress(payload)
    if "MS:1000523" in accessions:
        width, format_char = 8, "d"
    elif "MS:1000521" in accessions:
        width, format_char = 4, "f"
    else:
        raise ScoringError("mzML binary array lacks a supported 32- or 64-bit precision CV term.")
    if len(payload) % width:
        raise ScoringError("mzML binary array byte length is not divisible by its declared precision.")
    return kind, tuple(struct.unpack(f"<{len(payload) // width}{format_char}", payload))


def _peak_arrays(spectrum: ET.Element) -> tuple[tuple[float, ...], tuple[float, ...]]:
    mz: tuple[float, ...] = ()
    intensity: tuple[float, ...] = ()
    for element in spectrum.iter():
        if _local_name(element.tag) != "binaryDataArray":
            continue
        kind, values = _decode_array(element)
        if kind == "mz":
            mz = values
        elif kind == "intensity":
            intensity = values
    if len(mz) != len(intensity):
        raise ScoringError("mzML spectrum m/z and intensity arrays have different lengths.")
    if any(mz[index] > mz[index + 1] for index in range(len(mz) - 1)):
        raise ScoringError("mzML spectrum m/z array is not sorted.")
    return mz, intensity


def _nearest_index(mz: tuple[float, ...], target: float, ppm: float) -> int | None:
    if not mz:
        return None
    position = bisect.bisect_left(mz, target)
    candidates = [index for index in (position - 1, position) if 0 <= index < len(mz)]
    if not candidates:
        return None
    index = min(candidates, key=lambda candidate: abs(mz[candidate] - target))
    return index if abs(mz[index] - target) <= target * ppm * 1e-6 else None


def _envelope_indices(
    mz: tuple[float, ...], anchor_index: int, charge: int, ppm: float, max_isotopes: int
) -> tuple[int, ...]:
    spacing = ISOTOPE_MASS / charge
    anchor = mz[anchor_index]
    indices = {anchor_index}
    for direction in (-1, 1):
        for isotope in range(1, max_isotopes + 1):
            index = _nearest_index(mz, anchor + direction * isotope * spacing, ppm)
            if index is not None:
                indices.add(index)
    return tuple(sorted(indices))


def _competitor_envelopes(
    mz: tuple[float, ...], intensity: tuple[float, ...], window_indices: range,
    target_indices: set[int], ppm: float, max_isotopes: int
) -> list[tuple[tuple[int, ...], int, int]]:
    candidates: list[tuple[float, tuple[int, ...], int, int]] = []
    for index in window_indices:
        if index in target_indices or intensity[index] <= 0:
            continue
        for charge in range(1, 6):
            envelope = _envelope_indices(mz, index, charge, ppm, max_isotopes)
            envelope = tuple(candidate for candidate in envelope if candidate in window_indices and candidate not in target_indices)
            if len(envelope) >= 2:
                candidates.append((sum(intensity[candidate] for candidate in envelope), envelope, charge, index))
    selected: list[tuple[tuple[int, ...], int, int]] = []
    used: set[int] = set(target_indices)
    for _, envelope, charge, anchor_index in sorted(candidates, reverse=True):
        if any(index in used for index in envelope):
            continue
        selected.append((envelope, charge, anchor_index))
        used.update(envelope)
    return selected


def _prior_detection_count(scans: deque[Ms1Scan], target_mz: float, ppm: float) -> int:
    return sum(1 for scan in scans if _nearest_index(scan.mz, target_mz, ppm) is not None)


def _ms_level(params: dict[str, ET.Element]) -> int | None:
    value = _param(params, "MS:1000511", "ms level")
    return int(value.attrib["value"]) if value is not None else None


def _rt_seconds(params: dict[str, ET.Element]) -> float | None:
    parameter = _param(params, "MS:1000016", "scan start time")
    value = _float(parameter.attrib.get("value") if parameter is not None else None)
    if value is None or parameter is None:
        return value
    unit = (parameter.attrib.get("unitName") or parameter.attrib.get("unitAccession") or "").casefold()
    return value * 60 if unit in {"minute", "min", "uo:0000031"} else value


def _param_float(params: dict[str, ET.Element], accession: str, name: str) -> float | None:
    parameter = _param(params, accession, name)
    return _float(parameter.attrib.get("value") if parameter is not None else None)


def _classify(
    *, reason: str, interference: float | None, competitor_count: int,
    competitor_support: int | None, interference_threshold: float,
    minimum_competitor_prior_ms1_detections: int,
) -> tuple[str, str]:
    if reason:
        return "indeterminate", reason
    if interference is not None and interference >= interference_threshold and competitor_count > 0:
        if (competitor_support or 0) >= minimum_competitor_prior_ms1_detections:
            return "likely_chimeric", ""
        return "indeterminate", "competitor_insufficient_coelution_support"
    return "low_interference", ""


def score_bundle(
    bundle: RunBundle, *, interference_threshold: float, ppm_tolerance: float,
    max_isotopes: int, previous_ms1_scans: int,
    minimum_competitor_prior_ms1_detections: int,
) -> tuple[list[dict[str, str]], dict[str, str]]:
    """Score each mzML MS2 from RAW/mzML evidence only; PD data is not read."""
    if bundle.data_format != "mzml":
        raise ScoringError(f"Run {bundle.run_basename} requires mzML. Run `spacer convert` first.")
    prior_ms1: deque[Ms1Scan] = deque(maxlen=previous_ms1_scans)
    rows: list[dict[str, str]] = []
    counts = {"low_interference": 0, "likely_chimeric": 0, "indeterminate": 0}
    try:
        for _, spectrum in ET.iterparse(bundle.data_path, events=("end",)):
            if _local_name(spectrum.tag) != "spectrum":
                continue
            params = _spectrum_params(spectrum)
            level = _ms_level(params)
            native_id = spectrum.attrib.get("id", "")
            scan_id = _scan_number(native_id)
            if level == 1:
                mz, intensity = _peak_arrays(spectrum)
                prior_ms1.append(Ms1Scan(scan_id, _rt_seconds(params), mz, intensity))
                spectrum.clear()
                continue
            if level != 2:
                spectrum.clear()
                continue
            isolation_target = _param_float(params, "MS:1000827", "isolation window target m/z")
            lower_offset = _param_float(params, "MS:1000828", "isolation window lower offset")
            upper_offset = _param_float(params, "MS:1000829", "isolation window upper offset")
            charge_value = _param_float(params, "MS:1000041", "charge state")
            charge = int(charge_value) if charge_value and charge_value.is_integer() else None
            reason = ""
            target_fraction = interference = strongest_competitor = None
            competitor_count = 0
            target_support = competitor_support = None
            if not prior_ms1:
                reason = "no_preceding_ms1"
            elif isolation_target is None or lower_offset is None or upper_offset is None:
                reason = "missing_isolation_window_metadata"
            elif charge is None or charge <= 0:
                reason = "missing_or_invalid_precursor_charge"
            else:
                ms1 = prior_ms1[-1]
                window_start, window_end = isolation_target - lower_offset, isolation_target + upper_offset
                start = bisect.bisect_left(ms1.mz, window_start)
                end = bisect.bisect_right(ms1.mz, window_end)
                window_indices = range(start, end)
                total_intensity = sum(max(ms1.intensity[index], 0) for index in window_indices)
                target_index = _nearest_index(ms1.mz, isolation_target, ppm_tolerance)
                if target_index is None or target_index not in window_indices:
                    reason = "target_peak_not_found_in_preceding_ms1"
                elif total_intensity <= 0:
                    reason = "no_positive_intensity_in_isolation_window"
                else:
                    target_indices = set(_envelope_indices(ms1.mz, target_index, charge, ppm_tolerance, max_isotopes))
                    target_indices.intersection_update(window_indices)
                    target_intensity = sum(max(ms1.intensity[index], 0) for index in target_indices)
                    target_fraction = target_intensity / total_intensity
                    interference = 1 - target_fraction
                    competitors = _competitor_envelopes(ms1.mz, ms1.intensity, window_indices, target_indices, ppm_tolerance, max_isotopes)
                    competitor_count = len(competitors)
                    strongest_competitor = max((sum(ms1.intensity[index] for index in envelope) / total_intensity for envelope, _, _ in competitors), default=0.0)
                    target_support = _prior_detection_count(prior_ms1, ms1.mz[target_index], ppm_tolerance)
                    competitor_support = max((_prior_detection_count(prior_ms1, ms1.mz[anchor_index], ppm_tolerance) for _, _, anchor_index in competitors), default=0)
                    reason = ""
            classification, reason = _classify(
                reason=reason,
                interference=interference,
                competitor_count=competitor_count,
                competitor_support=competitor_support,
                interference_threshold=interference_threshold,
                minimum_competitor_prior_ms1_detections=minimum_competitor_prior_ms1_detections,
            )
            counts[classification] += 1
            rows.append({
                "run_basename": bundle.run_basename, "scan_id": str(scan_id or ""), "native_id": native_id,
                "ms2_rt_seconds": str(_rt_seconds(params) or ""),
                "preceding_ms1_scan_id": str(prior_ms1[-1].scan_id if prior_ms1 else ""),
                "precursor_mz": str(isolation_target or ""), "reported_charge": str(charge or ""),
                "isolation_lower_mz": str(isolation_target - lower_offset if isolation_target is not None and lower_offset is not None else ""),
                "isolation_upper_mz": str(isolation_target + upper_offset if isolation_target is not None and upper_offset is not None else ""),
                "target_precursor_fraction": str(target_fraction if target_fraction is not None else ""),
                "interference_fraction": str(interference if interference is not None else ""),
                "strongest_competitor_fraction": str(strongest_competitor if strongest_competitor is not None else ""),
                "competitor_count": str(competitor_count),
                "target_prior_ms1_detections": str(target_support) if target_support is not None else "",
                "competitor_prior_ms1_detections": str(competitor_support) if competitor_support is not None else "",
                "classification": classification, "indeterminate_reason": reason,
            })
            spectrum.clear()
    except (ET.ParseError, OSError, ValueError, zlib.error) as exc:
        raise ScoringError(f"Cannot score {bundle.data_path}: {exc}") from exc
    scorable = counts["low_interference"] + counts["likely_chimeric"]
    return rows, {
        "run_basename": bundle.run_basename, "total_ms2": str(len(rows)), "scorable_ms2": str(scorable),
        "low_interference": str(counts["low_interference"]), "likely_chimeric": str(counts["likely_chimeric"]),
        "indeterminate": str(counts["indeterminate"]),
        "likely_chimeric_fraction": str(counts["likely_chimeric"] / scorable if scorable else ""),
        "interference_threshold": str(interference_threshold), "ppm_tolerance": str(ppm_tolerance),
        "max_isotopes": str(max_isotopes), "previous_ms1_scans": str(previous_ms1_scans),
        "minimum_competitor_prior_ms1_detections": str(minimum_competitor_prior_ms1_detections),
    }


def sensitivity_summary(
    rows: list[dict[str, str]], *, thresholds: tuple[float, ...],
    minimum_competitor_prior_ms1_detections: int,
) -> list[dict[str, str]]:
    """Summarize the same independent evidence over predeclared thresholds."""
    summaries: list[dict[str, str]] = []
    for threshold in thresholds:
        counts = {"low_interference": 0, "likely_chimeric": 0, "indeterminate": 0}
        for row in rows:
            interference = _float(row["interference_fraction"])
            classification, _ = _classify(
                reason="" if interference is not None else row["indeterminate_reason"],
                interference=interference,
                competitor_count=int(row["competitor_count"]),
                competitor_support=int(row["competitor_prior_ms1_detections"] or 0),
                interference_threshold=threshold,
                minimum_competitor_prior_ms1_detections=minimum_competitor_prior_ms1_detections,
            )
            counts[classification] += 1
        scorable = counts["low_interference"] + counts["likely_chimeric"]
        summaries.append({
            "run_basename": rows[0]["run_basename"], "interference_threshold": str(threshold),
            "total_ms2": str(len(rows)), "scorable_ms2": str(scorable),
            "low_interference": str(counts["low_interference"]),
            "likely_chimeric": str(counts["likely_chimeric"]),
            "indeterminate": str(counts["indeterminate"]),
            "likely_chimeric_fraction": str(counts["likely_chimeric"] / scorable if scorable else ""),
            "minimum_competitor_prior_ms1_detections": str(minimum_competitor_prior_ms1_detections),
        })
    return summaries
