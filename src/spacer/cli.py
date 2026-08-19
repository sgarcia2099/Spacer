"""Command-line interface for Spacer's input-validation milestone."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path

from . import __version__
from .bundles import BundleValidationError, discover_bundles
from .conversion import ConversionError, convert_raw, existing_mzml_provenance, plan_conversions
from .reconciliation import ReconciliationError, reconcile_bundle


PWIZ_IMAGE = "proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses"


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="spacer",
        description="Validate Single-mode DDA input bundles for Spacer.",
    )
    parser.add_argument("--version", action="version", version=f"spacer {__version__}")
    commands = parser.add_subparsers(dest="command", required=True)

    doctor = commands.add_parser("doctor", help="Check the local foundation dependencies.")
    doctor.add_argument("--json", action="store_true", help="Emit machine-readable output.")

    for name, help_text in (
        ("single", "Validate one Single-mode bundle or one single-group directory."),
        ("validate", "Validate inputs and write a detailed input manifest."),
    ):
        command = commands.add_parser(name, help=help_text)
        source = command.add_mutually_exclusive_group(required=True)
        source.add_argument("--input", type=Path, help="One .raw or .mzML input file.")
        source.add_argument("--input-dir", type=Path, help="Flat directory of complete run bundles.")
        command.add_argument("--output-dir", type=Path, required=True, help="Directory for the input manifest.")
        command.add_argument(
            "--dry-run",
            action="store_true",
            help="Validate and print bundles without creating an output directory.",
        )
    convert = commands.add_parser(
        "convert", help="Convert only RAW files that lack a matching mzML file."
    )
    convert.add_argument("--input-dir", type=Path, required=True, help="Flat directory containing .raw/.mzML files.")
    convert.add_argument("--output-dir", type=Path, required=True, help="Directory for conversion_manifest.tsv.")
    convert.add_argument(
        "--dry-run",
        action="store_true",
        help="Print the conversion plan without creating mzML files or a manifest.",
    )
    reconcile = commands.add_parser(
        "reconcile", help="Reconcile converted mzML MS2 scans with Proteome Discoverer exports."
    )
    reconcile.add_argument("--input-dir", type=Path, required=True, help="Flat directory of converted Single-mode bundles.")
    reconcile.add_argument("--output-dir", type=Path, required=True, help="Directory for reconciliation TSV outputs.")
    reconcile.add_argument("--q-value-cutoff", type=float, default=0.01, help="PSM q-value cutoff used to mark an MS2 as identified (default: 0.01).")
    reconcile.add_argument("--rt-tolerance-seconds", type=float, default=1.0, help="Exact-scan metadata RT tolerance in seconds (default: 1.0).")
    reconcile.add_argument("--precursor-mz-tolerance-da", type=float, default=0.01, help="Exact-scan precursor m/z tolerance in Da (default: 0.01).")
    return parser


def _run_doctor(as_json: bool) -> int:
    docker_available = shutil.which("docker") is not None
    native_msconvert = shutil.which("msconvert") is not None
    container_msconvert = False
    container_detail = "image not checked"

    # Doctor must remain read-only: do not let `docker run` pull a missing image.
    if docker_available and not native_msconvert:
        try:
            image_check = subprocess.run(
                ["docker", "image", "inspect", PWIZ_IMAGE],
                capture_output=True,
                text=True,
                check=False,
                timeout=15,
            )
        except (OSError, subprocess.TimeoutExpired):
            image_check = None
            container_detail = "Docker image check failed or timed out"
        if image_check is not None and image_check.returncode == 0:
            try:
                converter_check = subprocess.run(
                    ["docker", "run", "--rm", PWIZ_IMAGE, "wine", "msconvert", "--help"],
                    capture_output=True,
                    text=True,
                    check=False,
                    timeout=60,
                )
                container_msconvert = converter_check.returncode == 0
                container_detail = (
                    "verified in local ProteoWizard container"
                    if container_msconvert
                    else "local ProteoWizard container could not run msconvert"
                )
            except (OSError, subprocess.TimeoutExpired):
                container_detail = "local ProteoWizard container check failed or timed out"
        elif image_check is not None:
            container_detail = "ProteoWizard container image is not installed locally"

    msconvert_available = native_msconvert or container_msconvert
    if native_msconvert:
        msconvert_source = "host PATH"
    elif container_msconvert:
        msconvert_source = "local ProteoWizard container"
    else:
        msconvert_source = container_detail

    checks = {
        "python": {"available": sys.version_info >= (3, 11), "version": sys.version.split()[0]},
        "docker": {"available": docker_available},
        "msconvert": {"available": msconvert_available, "source": msconvert_source},
    }
    if as_json:
        print(json.dumps(checks, indent=2, sort_keys=True))
    else:
        for name, check in checks.items():
            status = "ok" if check["available"] else "not found"
            version = f" ({check['version']})" if "version" in check else ""
            source = f" [{check['source']}]" if "source" in check else ""
            print(f"{name}: {status}{version}{source}")
        print("Docker/msconvert are required only when a supplied RAW needs conversion.")
    return 0 if checks["python"]["available"] else 2


def _write_manifest(output_dir: Path, rows: list[dict[str, str]]) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = output_dir / "run_manifest.tsv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return manifest


def _run_bundle_validation(args: argparse.Namespace) -> int:
    try:
        bundles = discover_bundles(input_path=args.input, input_dir=args.input_dir)
    except BundleValidationError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    rows = [bundle.manifest_row() for bundle in bundles]
    print(f"Validated {len(rows)} complete Single-mode run bundle(s):")
    for row in rows:
        pd_status = "available" if row["pd_isolation_interference_available"] == "true" else "not exported"
        print(f"- {row['run_basename']} ({row['data_format']}; PD interference: {pd_status})")
    if args.dry_run:
        print("Dry run: no files were created.")
    else:
        manifest = _write_manifest(args.output_dir, rows)
        print(f"Wrote input manifest: {manifest}")
    print("Input validation is complete. RAW conversion and chimericity scoring are not implemented in this milestone.")
    return 0


def _write_conversion_manifest(output_dir: Path, rows: list[dict[str, str]]) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = output_dir / "conversion_manifest.tsv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return manifest


def _run_conversion(args: argparse.Namespace) -> int:
    try:
        plans = plan_conversions(args.input_dir)
    except ConversionError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    pending = [plan for plan in plans if plan.status == "pending_conversion"]
    skipped = [plan for plan in plans if plan.status == "skipped_existing_mzml"]
    print(f"Conversion plan: {len(pending)} pending, {len(skipped)} existing mzML file(s) to validate.")
    for plan in plans:
        source = plan.raw_path.name if plan.raw_path is not None else "no RAW counterpart"
        print(f"- {plan.run_basename}: {plan.status} ({source} -> {plan.mzml_path.name})")
    if args.dry_run:
        print("Dry run: no Docker command, mzML validation, checksum, or file write was performed.")
        return 0

    rows: list[dict[str, str]] = []
    for plan in plans:
        try:
            if plan.status == "pending_conversion":
                summary, provenance = convert_raw(plan)
                status = "converted_and_validated"
            else:
                summary, provenance = existing_mzml_provenance(plan)
                status = "existing_mzml_validated"
        except ConversionError as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 2
        rows.append(
            {
                "run_basename": plan.run_basename,
                "status": status,
                "raw_path": str(plan.raw_path) if plan.raw_path is not None else "",
                "mzml_path": str(plan.mzml_path),
                "spectrum_count": str(summary.spectrum_count),
                "ms1_count": str(summary.ms1_count),
                "ms2_count": str(summary.ms2_count),
                **provenance,
            }
        )
    manifest = _write_conversion_manifest(args.output_dir, rows)
    print(f"Conversion and structural mzML validation complete. Wrote: {manifest}")
    return 0


def _write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _run_reconciliation(args: argparse.Namespace) -> int:
    if args.q_value_cutoff < 0:
        print("error: --q-value-cutoff must be non-negative.", file=sys.stderr)
        return 2
    if args.rt_tolerance_seconds < 0 or args.precursor_mz_tolerance_da < 0:
        print("error: reconciliation tolerances must be non-negative.", file=sys.stderr)
        return 2
    try:
        bundles = discover_bundles(input_path=None, input_dir=args.input_dir)
        all_rows: list[dict[str, str]] = []
        summaries: list[dict[str, str]] = []
        for bundle in bundles:
            rows, summary = reconcile_bundle(
                bundle,
                q_value_cutoff=args.q_value_cutoff,
                rt_tolerance_seconds=args.rt_tolerance_seconds,
                precursor_mz_tolerance_da=args.precursor_mz_tolerance_da,
            )
            all_rows.extend(rows)
            summaries.append(summary)
    except (BundleValidationError, ReconciliationError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    args.output_dir.mkdir(parents=True, exist_ok=True)
    _write_tsv(args.output_dir / "scan_reconciliation.tsv", all_rows)
    _write_tsv(args.output_dir / "scan_reconciliation_summary.tsv", summaries)
    exact_matches = sum(int(summary["matched_metadata_exact_precursor"]) for summary in summaries)
    isolation_window_matches = sum(
        int(summary["matched_metadata_precursor_within_isolation_window"])
        for summary in summaries
    )
    mismatched = sum(int(summary["matched_scan_metadata_mismatch_or_missing"]) for summary in summaries)
    print(
        f"Reconciled {len(all_rows)} scan keys across {len(summaries)} run(s): "
        f"{exact_matches} exact precursor matches, "
        f"{isolation_window_matches} precursor-in-window matches, "
        f"{mismatched} metadata mismatch/missing."
    )
    print(f"Wrote: {args.output_dir / 'scan_reconciliation.tsv'}")
    print(f"Wrote: {args.output_dir / 'scan_reconciliation_summary.tsv'}")
    print("Reconciliation is complete. Spacer chimericity scoring is not implemented in this milestone.")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    if args.command == "doctor":
        return _run_doctor(args.json)
    if args.command == "convert":
        return _run_conversion(args)
    if args.command == "reconcile":
        return _run_reconciliation(args)
    return _run_bundle_validation(args)


if __name__ == "__main__":
    raise SystemExit(main())
