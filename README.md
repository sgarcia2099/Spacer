# Spacer

Spacer is a bottom-up LC-MS/MS DDA chimericity workflow. The currently
available commands validate Single-mode input bundles, convert unmatched Thermo
RAW files to structurally validated mzML files, and reconcile converted MS2
metadata with Proteome Discoverer exports.

Project design, staged implementation notes, and planned analysis behavior are
in [misc/PLAN.md](misc/PLAN.md).

## Single-mode inputs

A Single-mode group contains one run or multiple technical replicates from the
same acquisition method. Each run needs one RAW or mzML input and three
tab-separated Proteome Discoverer exports sharing the same run stem:

```text
<run>.raw or <run>.mzML
<run>_MSMSSpectrumInfo.txt
<run>_PSMs.txt
<run>_PeptideGroups.txt
```

For example, `example_data/` contains three complete SC1 technical-replicate
bundles:

```text
02022026_SynComBench_QE1_180minGrad_100ulmin_5ug_SC1_R1.*
02022026_SynComBench_QE1_180minGrad_100ulmin_5ug_SC1_R2.*
02022026_SynComBench_QE1_180minGrad_100ulmin_5ug_SC1_R3.*
```

The input validator checks that each bundle has the three matching exports and
that their required headers are present:

| Export | Required headers |
| --- | --- |
| `*_MSMSSpectrumInfo.txt` | `First Scan`, `MS Order` |
| `*_PSMs.txt` | `First Scan`, `File ID`, `Annotated Sequence`, `q-Value`, `PEP` |
| `*_PeptideGroups.txt` | `Annotated Sequence`, `Number of PSMs` |

`Isolation Interference in Percent` is optional in the MS/MS spectrum-info
export. Its availability is recorded in the input manifest for the later,
separate agreement analysis.

## Install

```bash
git clone <REPOSITORY_URL> Spacer
cd Spacer
python3.11 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e .
```

Confirm the local foundation:

```bash
spacer doctor
```

`doctor` reports Python, Docker, and `msconvert` availability. Docker and
`msconvert` are only required when a RAW file later needs conversion to mzML.

## Convert unmatched RAW files

Review the conversion plan first. This command does not access Docker or write
files:

```bash
spacer convert \
  --input-dir example_data \
  --output-dir /tmp/spacer-conversion \
  --dry-run
```

To perform the conversion, omit `--dry-run`:

```bash
spacer convert \
  --input-dir example_data \
  --output-dir results/SC1_conversion
```

For each RAW without a matching case-insensitive mzML basename, Spacer runs
vendor centroiding through the local ProteoWizard Docker image and writes
`<run>.mzML` beside the RAW file. Existing matching mzML files are never
overwritten. The command writes `conversion_manifest.tsv`, including the
immutable container image digest, input/output SHA-256 values, and MS1/MS2
spectrum counts from structural mzML validation.

## Reconcile scan metadata

After conversion, reconcile every mzML MS2 scan with the corresponding
Proteome Discoverer spectrum and PSM exports:

```bash
spacer reconcile \
  --input-dir example_data \
  --output-dir /tmp/spacer-reconciliation
```

This writes `scan_reconciliation.tsv` and
`scan_reconciliation_summary.tsv`. The join key is the exact raw-data scan
number; Spacer never substitutes retention time, m/z, or row order for a
missing scan match. It records three matched-scan outcomes:

- `matched_metadata_exact_precursor`: retention time, precursor m/z, and charge
  agree within the declared tolerances.
- `matched_metadata_precursor_within_isolation_window`: retention time and
  charge agree, and Proteome Discoverer's precursor is inside the recorded mzML
  isolation window. This preserves a valid scan match when the search result
  reports a monoisotopic precursor rather than the isolation target.
- `matched_scan_metadata_mismatch_or_missing`: an exact scan match exists but
  its available metadata is incompatible or incomplete.

The tables are an audit of metadata agreement only. They do not calculate
chimericity or use Proteome Discoverer interference to modify Spacer results.
When chimericity scoring is added, raw/mzML evidence will remain the source of
the Spacer score; Proteome Discoverer fields will be reference context and a
separate post hoc correlation only, never ground truth or a one-to-one target.

## Validate inputs

Validate the complete SC1 technical-replicate group without writing files:

```bash
spacer single \
  --input-dir example_data \
  --output-dir results/SC1 \
  --dry-run
```

For one run only:

```bash
spacer single \
  --input example_data/02022026_SynComBench_QE1_180minGrad_100ulmin_5ug_SC1_R1.raw \
  --output-dir results/SC1_R1 \
  --dry-run
```

Omit `--dry-run` to write the verified input manifest:

```bash
spacer validate \
  --input-dir example_data \
  --output-dir results/SC1_validation
```

The command writes:

```text
results/SC1_validation/run_manifest.tsv
```

The manifest records the selected raw/mzML source for each run, paths to the
three Proteome Discoverer exports, and whether the Proteome Discoverer isolation
interference field is available.

## Validation rules

- Input directories are flat; subdirectories are not searched.
- RAW/mzML matching is case-insensitive by filename stem.
- If both a RAW and matching mzML are present, mzML is selected; no RAW
  conversion occurs.
- Duplicate RAW, duplicate mzML, duplicate result exports, incomplete bundles,
  unreadable export files, duplicate headers, and missing required headers are
  errors.
- `--dry-run` creates no output directory or files.

Input validation never changes RAW files. Conversion creates only missing mzML
files and never overwrites or removes RAW inputs.
