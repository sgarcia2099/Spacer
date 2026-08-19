# Spacer implementation plan

## Project status

**Input validation, idempotent RAW-to-mzML conversion, and scan-level mzML /
Proteome Discoverer reconciliation are implemented; manual reconciliation
verification and chimericity scoring remain pending.**

The scope is now one experimental group. The repository has three Thermo RAW technical replicates and their Proteome Discoverer exports in `example_data/`. The repository now provides a packaged `spacer` CLI with `doctor`, `single`, `validate`, `convert`, and `reconcile` commands. It validates input bundles, converts only unmatched RAW files through ProteoWizard Docker, writes conversion provenance, structurally validates mzML, and reconciles exact MS2 scan numbers with Proteome Discoverer exports. It does not yet calculate chimericity.

## Scope: Single mode

Single mode accepts one complete run bundle or a directory of complete bundles from a single experimental group. The group may contain technical replicates, but all runs use the same acquisition method. Single mode must not accept method labels, calculate between-group statistics, or present method comparisons.

Each run bundle contains:

```text
<run>.raw or <run>.mzML
<run>_MSMSSpectrumInfo.txt
<run>_PSMs.txt
<run>_PeptideGroups.txt
```

The current `SC1` group has three technical replicates:

| Run stem | RAW size | Results files |
| --- | ---: | --- |
| `02022026_SynComBench_QE1_180minGrad_100ulmin_5ug_SC1_R1` | 942,583,982 B | Spectrum info, PSMs, peptide groups |
| `02022026_SynComBench_QE1_180minGrad_100ulmin_5ug_SC1_R2` | 946,962,145 B | Spectrum info, PSMs, peptide groups |
| `02022026_SynComBench_QE1_180minGrad_100ulmin_5ug_SC1_R3` | 955,345,641 B | Spectrum info, PSMs, peptide groups |

Discovery must use the exact run stem, reject incomplete or ambiguous bundles, and never treat the three replicates as separate methods.

## Proteome Discoverer result import

Proteome Discoverer data provides external context only. It must not enter the Spacer precursor-interference formula or alter Spacer classifications.

**Reference-not-ground-truth rule:** Spacer must derive its chimericity score
and classification from raw/mzML evidence under its own declared assumptions.
Proteome Discoverer fields may provide scan-linked context and a separate,
post hoc agreement/correlation analysis, but must never train, recalibrate,
override, or serve as the one-to-one expected target for Spacer results.

### `*_MSMSSpectrumInfo.txt`

This is the primary external scan table because it has identified and unidentified MS2 spectra. The current exports contain `First Scan`, `File ID`, `RT in min`, `MS Order`, `Number of PSMs`, precursor m/z/charge, and `Isolation Interference in Percent`.

Import requirements:

- Parse tab-separated fields with a unique header.
- Require an integer `First Scan`, MS2 rows, and unique scan keys per run.
- Validate RT, m/z, and charge against mzML within documented tolerances.
- Parse optional `Isolation Interference in Percent` to the 0–100 range.
- Retain missing or invalid optional values as missing and report their counts.

### `*_PSMs.txt`

The PSM export uses `First Scan` and `File ID` to provide scan-linked identification context. Current fields include sequence, modifications, charge, m/z, `q-Value`, `PEP`, `PSM Ambiguity`, `Total Number of PSMs`, precursor abundance, and apex RT.

Group all PSM rows by exact `First Scan`. Define `identified` as at least one
PSM at or below a configurable q-value cutoff (default 0.01); retain the raw
PSM count and best q-value in the output. Use this context to summarize
identification confidence and ambiguity, cross-check precursor m/z/charge, and
choose manual inspection candidates. Do not use PSM scores or engine-specific
interference values in the independent Spacer score.

### `*_PeptideGroups.txt`

Peptide groups provide sequence, modifications, PSM count, protein ambiguity, abundance, confidence summaries, and apex RT. The current exports do not have a direct scan ID. Use them as peptide-level context and validate parsable exact PSM-to-peptide associations. Do not construct scan-level joins from RT or m/z alone; any non-exact association must be labelled contextual.

## Independent chimericity calculation

The primary calculation uses RAW/mzML only. For every eligible MS2:

1. Link to the preceding MS1 survey scan.
2. Read target precursor m/z and exact asymmetric isolation-window bounds.
3. Reconstruct the target peptide isotope envelope using reported or inferred charge and isotope spacing `1.00335 / z`.
4. Include only target envelope signal transmitted through the isolation window.
5. Detect non-target, peptide-like isotope envelopes in the same window.
6. Use neighboring MS1 scans to establish co-elution.
7. Record target precursor fraction, interference fraction, strongest competitor fraction, competitor count, co-elution support, assignment confidence, and MS1-to-MS2 time gap.
8. Classify the scan as `low_interference`, `likely_chimeric`, or `indeterminate`.

The default `likely_chimeric` rule is non-target interference of at least 0.30 and at least one co-eluting peptide-like competitor. This must be configurable, and outputs must include a prespecified sensitivity grid. `indeterminate` scans must never silently count as low interference.

## Proteome Discoverer agreement analysis

After independent Spacer scoring, perform an optional descriptive comparison only if `Isolation Interference in Percent` is available:

1. Join by raw-data run and exact `First Scan`/mzML scan number. Do not fall
   back to retention time, precursor m/z, or row order when an exact scan join
   fails; retain every unmatched record in the reconciliation ledger.
2. Retain both values unchanged: `100 * interference_fraction` and the exported Proteome Discoverer percentage.
3. Report matched, missing, and invalid-value counts.
4. Calculate Pearson and Spearman correlations on matched finite values.
5. Calculate median signed and median absolute differences in percentage points.
6. Repeat descriptively for all matched MS2, identified MS2, and unidentified MS2 where PSM linkage permits.
7. List the most discordant scan IDs for optional inspection.

This is an agreement result, not accuracy, calibration, or ground truth. It must not modify Spacer scores, classifications, or thresholds.

## Modes

### Single

Processes one run bundle or one directory of technical replicates. Outputs per-run results plus descriptive group statistics: number of runs, per-run likely-chimeric fractions, median, range, and indeterminate fractions. No hypothesis testing or method comparisons are allowed.

### Validation

Runs the same biological calculation as Single and adds bundle-discovery, conversion, mzML, scan/metadata join, peptide-group association, envelope, Proteome Discoverer agreement, provenance, and log diagnostics. It must report failures rather than silently repair them. It never creates plots by default.

### Inspect

An explicit, scan-specific action creates an MS1 isolation-window stick view, MS2 stick view, and coordinate export for the chosen scan. It follows review of `inspection_candidates.tsv`. Candidate selection must include low-interference, likely-chimeric, threshold-adjacent, indeterminate, and Spacer-versus-Proteome Discoverer discordant scans.

## Conversion and validation

- Match RAW/mzML case-insensitively by basename and convert only unmatched RAW.
- Preserve RAW inputs untouched.
- Use ProteoWizard `msconvert` vendor centroiding, with peak picking before transformations.
- Record converter command, immutable container digest, versions, checksums, sizes, timestamps, scan counts, and exit status.
- Validate mzML readability, scan order, finite/non-negative intensities, sorted m/z arrays, MS level, precursor fields, and isolation-window bounds.
- Report fatal errors, skipped scans, warnings, and informational diagnostics separately.

## Explicit exclusions

- No database search: do not run SAGE or any other search engine.
- No method comparison: this release supports one experimental group only.
- No FLASHDeconv: it is not appropriate for bottom-up DDA peptide co-isolation scoring.
- No msPurity runtime dependency: use a native peptide-specific method rather than a metabolomics-oriented feature workflow.

## Implementation phases

1. **Complete:** Create package metadata, CLI, and `spacer doctor`.
2. **Complete:** Implement complete single-run bundle discovery and structural header validation.
3. **Implemented; pending manual conversion check:** Implement idempotent RAW-to-mzML conversion, provenance, and basic MS1/MS2 structural validation.
4. **Implemented; pending manual reconciliation check:** Parse and validate all three Proteome Discoverer exports and exact scan joins.
5. Implement memory-efficient mzML reading and independent envelope scoring.
6. Implement Single and Validation reports, including technical-replicate group summaries.
7. Implement Proteome Discoverer agreement tables and discordant-scan list.
8. Implement explicit inspection and manual-review import.
9. Add controlled tests, then run the three SC1 technical replicates only after validation succeeds.

## Completion checklist

- [x] Complete bundles are discovered for R1, R2, and R3; required export headers are present.
- [x] Raw and external-result scan IDs are reconciled with explicit accounting.
- [ ] Spacer scoring does not consume vendor interference values.
- [ ] Proteome Discoverer correlation is optional, scan-matched, and separate.
- [ ] Group summaries are descriptive only.
- [ ] Validation is complete and plotting remains opt-in.
