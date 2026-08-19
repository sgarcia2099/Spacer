"""Idempotent Thermo RAW to mzML conversion and structural validation."""

from __future__ import annotations

import hashlib
import subprocess
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path


PWIZ_IMAGE = "proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses"


class ConversionError(RuntimeError):
    """Raised when conversion prerequisites or output validation fail."""


@dataclass(frozen=True)
class ConversionPlan:
    run_basename: str
    raw_path: Path | None
    mzml_path: Path
    status: str


@dataclass(frozen=True)
class MzmlSummary:
    spectrum_count: int
    ms1_count: int
    ms2_count: int


def _data_by_stem(directory: Path, suffix: str) -> dict[str, Path]:
    files = [path for path in directory.iterdir() if path.is_file() and path.suffix.casefold() == suffix]
    output: dict[str, Path] = {}
    duplicates: list[str] = []
    for path in files:
        key = path.stem.casefold()
        if key in output:
            duplicates.append(f"{output[key].name}, {path.name}")
        else:
            output[key] = path
    if duplicates:
        raise ConversionError(
            f"Ambiguous {suffix} file(s) for the same basename: {'; '.join(duplicates)}."
        )
    return output


def plan_conversions(directory: Path) -> list[ConversionPlan]:
    """Plan RAW conversion, selecting existing matching mzML files when present."""

    if not directory.is_dir():
        raise ConversionError(f"Input directory does not exist: {directory}")
    raws = _data_by_stem(directory, ".raw")
    mzmls = _data_by_stem(directory, ".mzml")
    if not raws and not mzmls:
        raise ConversionError(f"No .raw or .mzML files found directly in: {directory}")

    plans: list[ConversionPlan] = []
    for key in sorted(set(raws) | set(mzmls)):
        raw = raws.get(key)
        mzml = mzmls.get(key)
        if mzml is not None:
            plans.append(
                ConversionPlan(
                    run_basename=mzml.stem,
                    raw_path=raw,
                    mzml_path=mzml,
                    status="skipped_existing_mzml",
                )
            )
        else:
            assert raw is not None
            plans.append(
                ConversionPlan(
                    run_basename=raw.stem,
                    raw_path=raw,
                    mzml_path=raw.with_suffix(".mzML"),
                    status="pending_conversion",
                )
            )
    return plans


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def inspect_mzml(path: Path) -> MzmlSummary:
    """Check that mzML is readable and count spectra by declared MS level.

    This intentionally avoids decoding peak arrays; peak-level validation belongs
    to the subsequent scoring stage.
    """

    if not path.is_file() or path.stat().st_size == 0:
        raise ConversionError(f"mzML output is missing or empty: {path}")

    spectrum_count = ms1_count = ms2_count = 0
    try:
        for _, element in ET.iterparse(path, events=("end",)):
            if element.tag.rsplit("}", 1)[-1] != "spectrum":
                continue
            spectrum_count += 1
            level = None
            for child in element.iter():
                if child.tag.rsplit("}", 1)[-1] != "cvParam":
                    continue
                if child.attrib.get("accession") == "MS:1000511" or child.attrib.get("name") == "ms level":
                    level = child.attrib.get("value")
                    break
            if level == "1":
                ms1_count += 1
            elif level == "2":
                ms2_count += 1
            element.clear()
    except (ET.ParseError, OSError) as exc:
        raise ConversionError(f"mzML is not readable XML: {path}: {exc}") from exc

    if spectrum_count == 0:
        raise ConversionError(f"mzML contains no spectra: {path}")
    if ms1_count == 0 or ms2_count == 0:
        raise ConversionError(
            f"mzML must contain MS1 and MS2 spectra: {path} (MS1={ms1_count}, MS2={ms2_count})"
        )
    return MzmlSummary(spectrum_count=spectrum_count, ms1_count=ms1_count, ms2_count=ms2_count)


def ensure_container_converter() -> str:
    """Require a locally installed and runnable ProteoWizard image without pulling it."""

    try:
        image_check = subprocess.run(
            ["docker", "image", "inspect", PWIZ_IMAGE],
            capture_output=True,
            text=True,
            check=False,
            timeout=15,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise ConversionError(f"Could not check Docker image availability: {exc}") from exc
    if image_check.returncode != 0:
        raise ConversionError(
            f"ProteoWizard image is not installed locally: {PWIZ_IMAGE}. "
            f"Run: docker pull {PWIZ_IMAGE}"
        )
    try:
        converter_check = subprocess.run(
            ["docker", "run", "--rm", PWIZ_IMAGE, "wine", "msconvert", "--help"],
            capture_output=True,
            text=True,
            check=False,
            timeout=60,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise ConversionError(f"Could not run msconvert in Docker: {exc}") from exc
    if converter_check.returncode != 0:
        raise ConversionError("The local ProteoWizard container could not run msconvert.")
    try:
        digest_check = subprocess.run(
            ["docker", "image", "inspect", "--format", "{{index .RepoDigests 0}}", PWIZ_IMAGE],
            capture_output=True,
            text=True,
            check=False,
            timeout=15,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise ConversionError(f"Could not resolve the ProteoWizard image digest: {exc}") from exc
    if digest_check.returncode != 0 or not digest_check.stdout.strip():
        raise ConversionError("Could not resolve the immutable ProteoWizard image digest.")
    return digest_check.stdout.strip()


def convert_raw(plan: ConversionPlan) -> tuple[MzmlSummary, dict[str, str]]:
    """Convert one planned RAW file in place and validate the newly written mzML."""

    if plan.status != "pending_conversion" or plan.raw_path is None:
        raise ConversionError(f"Conversion is not pending for {plan.run_basename}.")
    if plan.mzml_path.exists():
        raise ConversionError(
            f"Refusing to overwrite existing mzML output: {plan.mzml_path}"
        )

    image_digest = ensure_container_converter()
    data_dir = plan.raw_path.parent.resolve()
    docker_command = [
        "docker",
        "run",
        "--rm",
        "-v",
        f"{data_dir}:/data",
        PWIZ_IMAGE,
        "wine",
        "msconvert",
        f"/data/{plan.raw_path.name}",
        "--mzML",
        "--zlib",
        "--filter",
        "peakPicking true 1-",
        "-o",
        "/data",
    ]
    try:
        completed = subprocess.run(
            docker_command,
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError as exc:
        raise ConversionError(f"Could not start RAW conversion for {plan.raw_path}: {exc}") from exc
    if completed.returncode != 0:
        stderr = completed.stderr.strip() or "no stderr captured"
        raise ConversionError(
            f"msconvert failed for {plan.raw_path.name} (exit {completed.returncode}): {stderr}"
        )
    summary = inspect_mzml(plan.mzml_path)
    provenance = {
        "converter_image": PWIZ_IMAGE,
        "converter_image_digest": image_digest,
        "converter_command": " ".join(docker_command),
        "raw_sha256": _sha256(plan.raw_path),
        "mzml_sha256": _sha256(plan.mzml_path),
    }
    return summary, provenance


def existing_mzml_provenance(plan: ConversionPlan) -> tuple[MzmlSummary, dict[str, str]]:
    """Validate an existing mzML and collect the read-only provenance available."""

    summary = inspect_mzml(plan.mzml_path)
    provenance = {
        "converter_image": "not_available_existing_mzml",
        "converter_image_digest": "not_available_existing_mzml",
        "converter_command": "not_available_existing_mzml",
        "raw_sha256": _sha256(plan.raw_path) if plan.raw_path is not None else "not_available",
        "mzml_sha256": _sha256(plan.mzml_path),
    }
    return summary, provenance
