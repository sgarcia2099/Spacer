"""Discovery and structural validation for Single-mode input bundles."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


RESULT_SUFFIXES = {
    "msms_spectrum_info": "_MSMSSpectrumInfo.txt",
    "psms": "_PSMs.txt",
    "peptide_groups": "_PeptideGroups.txt",
}

REQUIRED_HEADERS = {
    "msms_spectrum_info": {"First Scan", "MS Order"},
    "psms": {"First Scan", "File ID", "Annotated Sequence", "q-Value", "PEP"},
    "peptide_groups": {"Annotated Sequence", "Number of PSMs"},
}

OPTIONAL_HEADERS = {
    "msms_spectrum_info": {"Isolation Interference in Percent"},
}


class BundleValidationError(ValueError):
    """Raised when a Single-mode input directory is incomplete or ambiguous."""


@dataclass(frozen=True)
class RunBundle:
    """One raw/mzML input and its three Proteome Discoverer exports."""

    run_basename: str
    data_path: Path
    data_format: str
    msms_spectrum_info: Path
    psms: Path
    peptide_groups: Path
    has_pd_isolation_interference: bool

    def manifest_row(self) -> dict[str, str]:
        return {
            "run_basename": self.run_basename,
            "status": "validated_input",
            "data_path": str(self.data_path),
            "data_format": self.data_format,
            "msms_spectrum_info_path": str(self.msms_spectrum_info),
            "psms_path": str(self.psms),
            "peptide_groups_path": str(self.peptide_groups),
            "pd_isolation_interference_available": str(
                self.has_pd_isolation_interference
            ).lower(),
        }


def _casefold_stem(path: Path) -> str:
    return path.stem.casefold()


def _read_headers(path: Path) -> set[str]:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            headers = next(reader, None)
    except OSError as exc:
        raise BundleValidationError(f"Cannot read {path}: {exc}") from exc

    if not headers:
        raise BundleValidationError(f"{path} is empty; expected a tab-separated header.")

    cleaned = [header.strip() for header in headers]
    if any(not header for header in cleaned):
        raise BundleValidationError(f"{path} has an empty header name.")
    duplicates = sorted({header for header in cleaned if cleaned.count(header) > 1})
    if duplicates:
        raise BundleValidationError(
            f"{path} has duplicate header(s): {', '.join(duplicates)}."
        )
    return set(cleaned)


def _validate_result_file(kind: str, path: Path) -> bool:
    headers = _read_headers(path)
    missing = sorted(REQUIRED_HEADERS[kind] - headers)
    if missing:
        raise BundleValidationError(
            f"{path} is missing required column(s): {', '.join(missing)}."
        )
    return OPTIONAL_HEADERS.get(kind, set()).issubset(headers)


def _result_kind_and_stem(path: Path) -> tuple[str, str] | None:
    filename = path.name
    for kind, suffix in RESULT_SUFFIXES.items():
        if filename.casefold().endswith(suffix.casefold()):
            return kind, filename[: -len(suffix)]
    return None


def _one_per_key(paths: Iterable[Path], label: str) -> dict[str, Path]:
    output: dict[str, Path] = {}
    duplicates: list[str] = []
    for path in paths:
        key = _casefold_stem(path)
        if key in output:
            duplicates.append(f"{output[key].name}, {path.name}")
        else:
            output[key] = path
    if duplicates:
        raise BundleValidationError(
            f"Ambiguous {label} file(s) for the same basename: {'; '.join(duplicates)}."
        )
    return output


def _select_data_files(directory: Path) -> dict[str, Path]:
    raw_files = [path for path in directory.iterdir() if path.is_file() and path.suffix.casefold() == ".raw"]
    mzml_files = [path for path in directory.iterdir() if path.is_file() and path.suffix.casefold() == ".mzml"]
    raws = _one_per_key(raw_files, "RAW")
    mzmls = _one_per_key(mzml_files, "mzML")

    selected = dict(raws)
    # A pre-existing mzML is authoritative and prevents RAW conversion later.
    selected.update(mzmls)
    return selected


def _discover_from_directory(directory: Path) -> list[RunBundle]:
    if not directory.is_dir():
        raise BundleValidationError(f"Input directory does not exist: {directory}")

    data_files = _select_data_files(directory)
    if not data_files:
        raise BundleValidationError(
            f"No .raw or .mzML files found directly in input directory: {directory}"
        )

    result_maps: dict[str, dict[str, Path]] = {kind: {} for kind in RESULT_SUFFIXES}
    duplicate_results: list[str] = []
    for path in directory.iterdir():
        if not path.is_file():
            continue
        parsed = _result_kind_and_stem(path)
        if parsed is None:
            continue
        kind, stem = parsed
        key = stem.casefold()
        if key in result_maps[kind]:
            duplicate_results.append(
                f"{result_maps[kind][key].name}, {path.name} ({kind})"
            )
        else:
            result_maps[kind][key] = path
    if duplicate_results:
        raise BundleValidationError(
            "Ambiguous Proteome Discoverer export(s): " + "; ".join(duplicate_results)
        )

    bundles: list[RunBundle] = []
    errors: list[str] = []
    for key, data_path in sorted(data_files.items()):
        missing = [kind for kind, paths in result_maps.items() if key not in paths]
        if missing:
            errors.append(
                f"{data_path.name}: missing result export(s): {', '.join(missing)}"
            )
            continue
        try:
            has_pd_interference = _validate_result_file(
                "msms_spectrum_info", result_maps["msms_spectrum_info"][key]
            )
            _validate_result_file("psms", result_maps["psms"][key])
            _validate_result_file("peptide_groups", result_maps["peptide_groups"][key])
        except BundleValidationError as exc:
            errors.append(str(exc))
            continue
        bundles.append(
            RunBundle(
                run_basename=data_path.stem,
                data_path=data_path,
                data_format=data_path.suffix.casefold().lstrip("."),
                msms_spectrum_info=result_maps["msms_spectrum_info"][key],
                psms=result_maps["psms"][key],
                peptide_groups=result_maps["peptide_groups"][key],
                has_pd_isolation_interference=has_pd_interference,
            )
        )
    if errors:
        raise BundleValidationError("Input bundle validation failed:\n- " + "\n- ".join(errors))
    return bundles


def discover_bundles(*, input_path: Path | None, input_dir: Path | None) -> list[RunBundle]:
    """Discover and structurally validate one or more Single-mode run bundles."""

    if (input_path is None) == (input_dir is None):
        raise BundleValidationError("Specify exactly one of --input or --input-dir.")
    if input_path is not None:
        if not input_path.is_file():
            raise BundleValidationError(f"Input file does not exist: {input_path}")
        if input_path.suffix.casefold() not in {".raw", ".mzml"}:
            raise BundleValidationError(
                f"Input must end in .raw or .mzML, not: {input_path.name}"
            )
        bundles = _discover_from_directory(input_path.parent)
        matching = [
            bundle
            for bundle in bundles
            if bundle.run_basename.casefold() == input_path.stem.casefold()
        ]
        if not matching:
            raise BundleValidationError(
                f"No complete result bundle found for input file: {input_path.name}"
            )
        return matching
    return _discover_from_directory(input_dir)
