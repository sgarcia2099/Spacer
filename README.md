# Spacer

Spacer is a bottom-up LC-MS/MS DDA chimericity workflow. The currently
available commands validate Single-mode input bundles, convert unmatched Thermo
RAW files to structurally validated mzML files, reconcile converted MS2
metadata with Proteome Discoverer exports, calculate independent MS1-based
precursor-interference scores, create descriptive Analysis and Validation
reports, and optionally inspect or compare exact-scan PD reference context.

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
python3 -m venv .venv
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
Raw/mzML evidence is the source of the Spacer score; Proteome Discoverer fields
remain reference context and will support a separate post hoc correlation only,
never ground truth or a one-to-one target.

## Score independent precursor interference

Score the converted mzML files from their MS1 peak arrays and MS2 isolation
metadata. This command does not read Proteome Discoverer values for the score:

```bash
spacer score \
  --input-dir example_data \
  --output-dir /tmp/spacer-scoring \
  --interference-threshold 0.30 \
  --ppm-tolerance 15 \
  --max-isotopes 6 \
  --previous-ms1-scans 3 \
  --minimum-competitor-prior-ms1-detections 2 \
  --sensitivity-thresholds 0.20,0.30,0.40
```

It writes `ms2_chimericity.tsv`, `run_summary.tsv`, and
`sensitivity_summary.tsv`. Each MS2 is labeled
`low_interference`, `likely_chimeric`, or `indeterminate`. The current scorer
uses only a preceding MS1 survey scan, the exact recorded isolation window,
reported precursor charge, and charge-aware isotope envelopes. It retains
prior-MS1 detection counts. By default, a competing envelope must be detected
in at least two of the three retained MS1 scans to be called co-eluting; an
otherwise high-interference scan without that support is `indeterminate`, not
low interference. The sensitivity table reports the same evidence at the
prespecified 0.20, 0.30, and 0.40 thresholds. It does not use Proteome
Discoverer scores or labels.

## Inspect selected spectra

Inspection is opt-in. Spacer never generates plots during scoring or candidate
selection. First create a compact, auditable candidate table from a completed
score:

~~~bash
spacer inspect \
  --input-dir example_data \
  --scoring-table /tmp/spacer-scoring-coelution/ms2_chimericity.tsv \
  --output-dir /tmp/spacer-inspection \
  --select-candidates \
  --balance-runs \
  --per-category 5
~~~

This writes inspection_candidates.tsv with five representatives, where
available, from each of: high interference, low interference,
threshold-adjacent, and indeterminate MS2. The balance-runs option ensures
each technical replicate contributes before any run is repeated. Review the
table before choosing an explicit scan:

~~~bash
column -ts $'\t' /tmp/spacer-inspection/inspection_candidates.tsv
~~~

Then export coordinates for one selected run/scan pair. Omit the plot option
to export coordinates only:

~~~bash
spacer inspect \
  --input-dir example_data \
  --scoring-table /tmp/spacer-scoring-coelution/ms2_chimericity.tsv \
  --output-dir /tmp/spacer-inspection \
  --run-basename 02022026_SynComBench_QE1_180minGrad_100ulmin_5ug_SC1_R2 \
  --scan-id 76074
~~~

The command writes an MS1 isolation-window coordinate table (including
isolation-window and target-envelope flags), the complete MS2 coordinate
table, and inspection metadata. Use those TSVs in any plotting environment of
your choice to verify that low-interference cases look clean and likely
chimeric cases visibly contain competing precursor signal.

To render an explicit PNG and SVG after choosing a scan, install the optional
plotting extra once:

~~~bash
python -m pip install -e '.[plot]'
~~~

Then add --plot to the scan-specific command. The plot distinguishes target
envelope signal inside the nominal isolation window from target-envelope peaks
outside it. Only target signal inside the amber isolation-window region enters
the current strict scoring numerator.

## Analysis mode

Create the concise, descriptive report for one technical-replicate group after
scoring. It reports per-run results, pooled counts, replicate median/range, and
the prespecified sensitivity grid. It performs no method comparison or
hypothesis test.

~~~bash
spacer analyze \
  --scoring-dir /tmp/spacer-scoring-coelution \
  --output-dir results/SC1_analysis
~~~

Outputs are analysis_summary.tsv, sensitivity_group_summary.tsv, and
analysis_report.md.

## Validation mode

Validation mode performs the same input/result accounting with detailed
diagnostics. It checks bundle discovery, mzML MS1/MS2 counts, unique score
scan IDs, score arithmetic, finite score values, sensitivity-grid coverage,
and optional reconciliation/agreement totals. It never recalibrates a Spacer
score or creates a plot.

~~~bash
spacer validation \
  --input-dir example_data \
  --scoring-dir /tmp/spacer-scoring-coelution \
  --reconciliation-dir /tmp/spacer-reconciliation \
  --agreement-dir /tmp/spacer-agreement \
  --output-dir results/SC1_validation
~~~

It writes validation_checks.tsv and validation_report.md.

## Optional Proteome Discoverer agreement

This step is separate from scoring. It joins only exact run/scan matches,
retains both interference values unchanged, and reports descriptive agreement
for all matched MS2 plus identified and unidentified subsets. It does not
train, alter, or validate the Spacer classifications as ground truth.

~~~bash
spacer agreement \
  --input-dir example_data \
  --scoring-dir /tmp/spacer-scoring-coelution \
  --output-dir results/SC1_pd_agreement \
  --q-value-cutoff 0.01 \
  --top-discordant 20
~~~

Outputs are pd_agreement.tsv, pd_agreement_summary.tsv, and
pd_discordant_candidates.tsv. Pass the first file to spacer inspect with
--agreement-table to add PD-discordant scans to an opt-in candidate table.

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
