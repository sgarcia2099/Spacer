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
from .agreement import AgreementError, agreement_for_bundles
from .bundles import BundleValidationError, discover_bundles
from .conversion import ConversionError, convert_raw, existing_mzml_provenance, plan_conversions
from .inspection import (
    InspectionError,
    export_scan_coordinates,
    read_discordant_rows,
    read_score_rows,
    render_scan_plot,
    select_candidates,
)
from .reconciliation import ReconciliationError, reconcile_bundle
from .reporting import ReportingError, build_analysis, write_tsv
from .scoring import ScoringError, score_bundle, sensitivity_summary
from .validation import ValidationError, validate_completed_analysis, write_validation


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
    score = commands.add_parser("score", help="Calculate independent MS1-based precursor interference scores.")
    score.add_argument("--input-dir", type=Path, required=True, help="Flat directory of converted Single-mode bundles.")
    score.add_argument("--output-dir", type=Path, required=True, help="Directory for independent scoring TSV outputs.")
    score.add_argument("--interference-threshold", type=float, default=0.30, help="Likely-chimeric interference threshold (default: 0.30).")
    score.add_argument("--ppm-tolerance", type=float, default=15.0, help="MS1 peak matching tolerance in ppm (default: 15).")
    score.add_argument("--max-isotopes", type=int, default=6, help="Maximum isotope offsets examined per envelope (default: 6).")
    score.add_argument("--previous-ms1-scans", type=int, default=3, help="Previous MS1 scans retained for persistence evidence (default: 3).")
    score.add_argument("--minimum-competitor-prior-ms1-detections", type=int, default=2, help="MS1 detections required for a competing envelope to be considered co-eluting (default: 2).")
    score.add_argument("--sensitivity-thresholds", default="0.20,0.30,0.40", help="Comma-separated, prespecified interference thresholds for sensitivity_summary.tsv (default: 0.20,0.30,0.40).")
    inspect = commands.add_parser(
        "inspect", help="Select manual-review candidates or export coordinates for one chosen MS2."
    )
    inspect.add_argument("--input-dir", type=Path, required=True, help="Flat directory of converted Single-mode bundles.")
    inspect.add_argument("--scoring-table", type=Path, required=True, help="ms2_chimericity.tsv produced by spacer score.")
    inspect.add_argument("--output-dir", type=Path, required=True, help="Directory for opt-in inspection outputs.")
    inspect.add_argument("--select-candidates", action="store_true", help="Write representative candidate rows; does not generate plots.")
    inspect.add_argument("--balance-runs", action="store_true", help="When selecting candidates, represent each run before selecting repeats.")
    inspect.add_argument("--agreement-table", type=Path, help="Optional pd_agreement.tsv used to add PD-discordant candidates.")
    inspect.add_argument("--plot", action="store_true", help="Render PNG and SVG for one explicitly selected scan.")
    inspect.add_argument("--run-basename", help="Run basename for one explicit coordinate export.")
    inspect.add_argument("--scan-id", type=int, help="MS2 scan number for one explicit coordinate export.")
    inspect.add_argument("--per-category", type=int, default=5, help="Candidate rows per category (default: 5).")
    inspect.add_argument("--interference-threshold", type=float, default=0.30, help="Score threshold used for threshold-adjacent candidates (default: 0.30).")
    inspect.add_argument("--ppm-tolerance", type=float, default=15.0, help="Envelope annotation tolerance in ppm (default: 15).")
    inspect.add_argument("--max-isotopes", type=int, default=6, help="Maximum isotope offsets for target annotation (default: 6).")
    inspect.add_argument("--window-margin-mz", type=float, default=0.5, help="Additional m/z shown on either side of the MS1 isolation window (default: 0.5).")
    analyze = commands.add_parser(
        "analyze", help="Create concise descriptive Single-mode reports from completed scoring outputs."
    )
    analyze.add_argument("--scoring-dir", type=Path, required=True, help="Directory containing score TSV outputs.")
    analyze.add_argument("--output-dir", type=Path, required=True, help="Directory for concise Analysis reports.")
    validation = commands.add_parser(
        "validation", help="Run detailed structural and arithmetic diagnostics on completed Single-mode outputs."
    )
    validation.add_argument("--input-dir", type=Path, required=True, help="Flat directory of converted Single-mode bundles.")
    validation.add_argument("--scoring-dir", type=Path, required=True, help="Directory containing score TSV outputs.")
    validation.add_argument("--reconciliation-dir", type=Path, help="Optional directory containing reconciliation TSV outputs.")
    validation.add_argument("--agreement-dir", type=Path, help="Optional directory containing PD agreement TSV outputs.")
    validation.add_argument("--output-dir", type=Path, required=True, help="Directory for Validation diagnostics.")
    agreement = commands.add_parser(
        "agreement", help="Describe exact-scan agreement with optional PD interference exports."
    )
    agreement.add_argument("--input-dir", type=Path, required=True, help="Flat directory of converted Single-mode bundles.")
    agreement.add_argument("--scoring-dir", type=Path, required=True, help="Directory containing ms2_chimericity.tsv.")
    agreement.add_argument("--output-dir", type=Path, required=True, help="Directory for descriptive PD agreement outputs.")
    agreement.add_argument("--q-value-cutoff", type=float, default=0.01, help="PD PSM q-value cutoff for identified subsets (default: 0.01).")
    agreement.add_argument("--top-discordant", type=int, default=20, help="Discordant finite scans retained per run (default: 20).")
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
    print("Input validation is complete. Continue with `spacer convert`, `spacer reconcile`, and `spacer score` as applicable.")
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
    print("Reconciliation is complete. Use `spacer score` for independent mzML-based scoring.")
    return 0


def _run_scoring(args: argparse.Namespace) -> int:
    if not 0 <= args.interference_threshold <= 1:
        print("error: --interference-threshold must be between 0 and 1.", file=sys.stderr)
        return 2
    if args.ppm_tolerance <= 0 or args.max_isotopes < 1 or args.previous_ms1_scans < 1 or args.minimum_competitor_prior_ms1_detections < 1:
        print("error: ppm tolerance, max isotopes, previous MS1 scans, and minimum competitor detections must be positive.", file=sys.stderr)
        return 2
    try:
        thresholds = tuple(sorted({float(value.strip()) for value in args.sensitivity_thresholds.split(",")}))
        if not thresholds or any(not 0 <= threshold <= 1 for threshold in thresholds):
            raise ValueError
        if args.interference_threshold not in thresholds:
            print("error: --sensitivity-thresholds must include --interference-threshold.", file=sys.stderr)
            return 2
        bundles = discover_bundles(input_path=None, input_dir=args.input_dir)
        all_rows: list[dict[str, str]] = []
        summaries: list[dict[str, str]] = []
        sensitivity_rows: list[dict[str, str]] = []
        for bundle in bundles:
            print(f"Scoring {bundle.run_basename} from {bundle.data_path.name}...", flush=True)
            rows, summary = score_bundle(
                bundle,
                interference_threshold=args.interference_threshold,
                ppm_tolerance=args.ppm_tolerance,
                max_isotopes=args.max_isotopes,
                previous_ms1_scans=args.previous_ms1_scans,
                minimum_competitor_prior_ms1_detections=args.minimum_competitor_prior_ms1_detections,
            )
            all_rows.extend(rows)
            summaries.append(summary)
            sensitivity_rows.extend(sensitivity_summary(
                rows,
                thresholds=thresholds,
                minimum_competitor_prior_ms1_detections=args.minimum_competitor_prior_ms1_detections,
            ))
    except (BundleValidationError, ScoringError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    args.output_dir.mkdir(parents=True, exist_ok=True)
    _write_tsv(args.output_dir / "ms2_chimericity.tsv", all_rows)
    _write_tsv(args.output_dir / "run_summary.tsv", summaries)
    _write_tsv(args.output_dir / "sensitivity_summary.tsv", sensitivity_rows)
    print(f"Independently scored {len(all_rows)} MS2 scans across {len(summaries)} run(s).")
    print(f"Wrote: {args.output_dir / 'ms2_chimericity.tsv'}")
    print(f"Wrote: {args.output_dir / 'run_summary.tsv'}")
    print(f"Wrote: {args.output_dir / 'sensitivity_summary.tsv'}")
    print("Proteome Discoverer values were not read for this scoring calculation.")
    return 0


def _run_inspection(args: argparse.Namespace) -> int:
    if not 0 <= args.interference_threshold <= 1:
        print("error: --interference-threshold must be between 0 and 1.", file=sys.stderr)
        return 2
    if args.per_category < 1 or args.ppm_tolerance <= 0 or args.max_isotopes < 1 or args.window_margin_mz < 0:
        print("error: inspection counts/tolerances must be positive (window margin may be zero).", file=sys.stderr)
        return 2
    choosing_candidates = args.select_candidates
    choosing_scan = args.run_basename is not None or args.scan_id is not None
    if choosing_candidates == choosing_scan:
        print("error: select either --select-candidates or both --run-basename and --scan-id.", file=sys.stderr)
        return 2
    if choosing_scan and (args.run_basename is None or args.scan_id is None):
        print("error: coordinate export requires both --run-basename and --scan-id.", file=sys.stderr)
        return 2
    try:
        rows = read_score_rows(args.scoring_table)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        if choosing_candidates:
            if args.plot:
                print("error: --plot requires one explicit --run-basename/--scan-id pair.", file=sys.stderr)
                return 2
            candidates = select_candidates(
                rows,
                interference_threshold=args.interference_threshold,
                per_category=args.per_category,
                balance_runs=args.balance_runs,
                discordant_rows=read_discordant_rows(args.agreement_table) if args.agreement_table else None,
            )
            _write_tsv(args.output_dir / "inspection_candidates.tsv", candidates)
            print(f"Wrote {len(candidates)} opt-in inspection candidates: {args.output_dir / 'inspection_candidates.tsv'}")
            print("No plots were generated.")
            return 0
        matching_rows = [
            row for row in rows
            if row["run_basename"] == args.run_basename and row["scan_id"] == str(args.scan_id)
        ]
        if len(matching_rows) != 1:
            raise InspectionError(
                f"Expected one scoring row for run {args.run_basename!r}, scan {args.scan_id}; found {len(matching_rows)}."
            )
        bundles = discover_bundles(input_path=None, input_dir=args.input_dir)
        matching_bundles = [bundle for bundle in bundles if bundle.run_basename == args.run_basename]
        if len(matching_bundles) != 1:
            raise InspectionError(f"Could not find one converted mzML bundle for run {args.run_basename!r}.")
        paths = export_scan_coordinates(
            matching_bundles[0],
            score_row=matching_rows[0],
            output_dir=args.output_dir,
            ppm_tolerance=args.ppm_tolerance,
            max_isotopes=args.max_isotopes,
            window_margin_mz=args.window_margin_mz,
        )
        plot_paths = render_scan_plot(paths) if args.plot else None
    except (BundleValidationError, InspectionError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(f"Wrote MS1-window coordinates: {paths['ms1']}")
    print(f"Wrote MS2 coordinates: {paths['ms2']}")
    print(f"Wrote inspection metadata: {paths['metadata']}")
    if plot_paths is None:
        print("No plot was generated; rerun with --plot to render an explicit PNG and SVG.")
    else:
        print(f"Wrote inspection PNG: {plot_paths['png']}")
        print(f"Wrote inspection SVG: {plot_paths['svg']}")
    return 0


def _run_analysis(args: argparse.Namespace) -> int:
    try:
        analysis_rows, sensitivity_rows, report = build_analysis(args.scoring_dir)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        write_tsv(args.output_dir / "analysis_summary.tsv", analysis_rows)
        write_tsv(args.output_dir / "sensitivity_group_summary.tsv", sensitivity_rows)
        (args.output_dir / "analysis_report.md").write_text(report, encoding="utf-8")
    except ReportingError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(f"Wrote: {args.output_dir / 'analysis_summary.tsv'}")
    print(f"Wrote: {args.output_dir / 'sensitivity_group_summary.tsv'}")
    print(f"Wrote: {args.output_dir / 'analysis_report.md'}")
    print("Analysis is descriptive for one technical-replicate group; no method comparison was performed.")
    return 0


def _run_validation_mode(args: argparse.Namespace) -> int:
    try:
        bundles = discover_bundles(input_path=None, input_dir=args.input_dir)
        checks, report = validate_completed_analysis(
            bundles,
            args.scoring_dir,
            args.reconciliation_dir,
            args.agreement_dir,
        )
        write_validation(args.output_dir, checks, report)
    except (BundleValidationError, ValidationError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(f"Validation passed: {len(checks)} checks.")
    print(f"Wrote: {args.output_dir / 'validation_checks.tsv'}")
    print(f"Wrote: {args.output_dir / 'validation_report.md'}")
    print("Validation reports diagnostics only; it does not change Spacer scores or classifications.")
    return 0


def _run_agreement(args: argparse.Namespace) -> int:
    try:
        bundles = discover_bundles(input_path=None, input_dir=args.input_dir)
        rows, summaries, discordant = agreement_for_bundles(
            bundles,
            args.scoring_dir / "ms2_chimericity.tsv",
            q_value_cutoff=args.q_value_cutoff,
            top_discordant=args.top_discordant,
        )
        args.output_dir.mkdir(parents=True, exist_ok=True)
        write_tsv(args.output_dir / "pd_agreement.tsv", rows)
        write_tsv(args.output_dir / "pd_agreement_summary.tsv", summaries)
        write_tsv(args.output_dir / "pd_discordant_candidates.tsv", discordant)
    except (BundleValidationError, AgreementError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(f"Wrote: {args.output_dir / 'pd_agreement.tsv'}")
    print(f"Wrote: {args.output_dir / 'pd_agreement_summary.tsv'}")
    print(f"Wrote: {args.output_dir / 'pd_discordant_candidates.tsv'}")
    print("PD agreement is descriptive reference context only; it did not modify Spacer scores or classifications.")
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
    if args.command == "score":
        return _run_scoring(args)
    if args.command == "inspect":
        return _run_inspection(args)
    if args.command == "analyze":
        return _run_analysis(args)
    if args.command == "validation":
        return _run_validation_mode(args)
    if args.command == "agreement":
        return _run_agreement(args)
    return _run_bundle_validation(args)


if __name__ == "__main__":
    raise SystemExit(main())
