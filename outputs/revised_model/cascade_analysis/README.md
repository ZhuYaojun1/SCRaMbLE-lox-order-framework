# Rearrangement-cascade analysis data

This directory contains the audited outputs of deterministic state replay across all 88 formal-grid scenarios. The analysis tests whether inversions or duplications preceding a terminal Essentiality Gate-failing deletion changed the accessibility, proposal probability, or Gate consequence of that later deletion path.

## Release snapshot

- 1,596,885 accepted events inspected.
- 438,708 Gate-failing trajectories.
- 196,246 Gate-failing trajectories with at least one prior inversion or duplication (44.73%).
- 328,098 replayed non-deletion precursor events: 241,208 inversions and 86,890 duplications.
- 8,785 deterministically sampled precursor states evaluated by exact all-pair Gate-failing deletion hazard.
- 88 of 88 scenarios reconciled, with zero replay mismatches.
- 13 of 13 automated validation checks passed.

## Core data

`cascade_precursor_events.csv.xz` is the event-level source table. Each row represents an inversion or duplication before a later terminal Gate-failing deletion. It records the exact terminal-pair state before and after the precursor and, for the scenario-stratified sample, the total one-step Gate-failing deletion hazard across all active pairs.

`cascade_trajectory_effects.csv.xz` is the trajectory-level source table for Gate-failing histories containing at least one non-deletion precursor. It supports history composition, trajectory-level effect fractions, and confidence intervals.

The remaining CSV files are compact pooled, scenario-level, chromosome-level, timing, gene-association, threshold, and ranked-event summaries. `cascade_replay_audit.csv` provides scenario-by-scenario accounting closure. `cascade_analysis_config.json` records the replay configuration.

To repeat the state replay after regenerating the formal scenario event logs, run:

```bash
python scripts/revised_model/analyze_rearrangement_cascade.py \
  --event-root /path/to/formal/scenarios \
  --inputs-dir outputs/revised_model/inputs \
  --output-dir outputs/revised_model/cascade_analysis \
  --gate-name strict \
  --total-hazard-sample-per-type 50 \
  --seed 20260718

python scripts/revised_model/finalize_rearrangement_cascade.py \
  --output-dir outputs/revised_model/cascade_analysis \
  --inputs-dir outputs/revised_model/inputs \
  --seed 20260718
```

## Reports and manuscript assets

- `CASCADE_ANALYSIS_REPORT.md` gives the principal numerical findings and interpretation boundary.
- `CASCADE_ANALYSIS_VALIDATION.md` records all automated checks.
- `manuscript_assets/figures/` contains 400-dpi PNG plus PDF and SVG versions of Figures 10 and 11.
- `manuscript_assets/source_data/` contains the exact plotted source data and Table 11 source table.
- `manuscript_assets/SCRaMbLE_section_3_8_Figure_Table_package.docx` is a standalone placement package, not the full manuscript.

Run `python scripts/revised_model/build_cascade_manuscript_assets.py` to regenerate the figure/table package into `manuscript_assets/rebuilt/` without replacing the audited exports.

## Interpretation boundary

Deletion is the only event type that can immediately produce failure under the copy-number Essentiality Gate. The retrospective terminal-pair analysis is conditioned on trajectories that later failed the Gate and therefore quantifies path modification, not a population causal effect. The sampled all-pair analysis measures a one-step model hazard under the specified pair-sampling and event rules. Gate-failing trajectories are model states and are not experimentally observed dead-cell genomes.

## Integrity

`SHA256SUMS.csv` records the path, byte size, and SHA-256 digest of every file in this directory except the checksum file itself. The repository-wide `MANIFEST_SHA256.csv` provides the corresponding release-level inventory.
