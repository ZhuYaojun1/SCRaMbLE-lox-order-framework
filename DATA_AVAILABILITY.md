# Data availability and provenance

## Public release contents

The repository contains model code, model-ready derived inputs, compact main-grid source data, statistical outputs, robustness summaries, quality-control reports, and manuscript figures. These files support inspection of the numerical claims without redistributing raw sequencing data or large event-level logs.

| Evidence layer | Repository location | Purpose |
|---|---|---|
| Audited lox inputs | `outputs/revised_model/inputs/` | Active-site filtering, coordinates, Hi-C coverage, exact gene mapping, and Essentiality Gate membership |
| Main-grid source data | `outputs/revised_model/main_run/` | Scenario totals, endpoint counts, and site-level predicted frequencies |
| Statistical analyses | `outputs/revised_model/analysis/` | Rarefaction, Pareto frontier, backtesting, hotspot metrics, null controls, holdout, gene risk, and endpoint-alias audits |
| Robustness analyses | `outputs/revised_model/robustness/` | Five-seed and Essentiality Gate sensitivity summaries |
| Rearrangement-cascade analysis | `outputs/revised_model/cascade_analysis/` | Deterministic replay of precursor inversion/duplication histories, exact terminal-pair effects, sampled all-pair Gate-failing deletion hazard, validation, and Section 3.8 manuscript assets |
| Quality control | `outputs/revised_model/reports/` | Dry-run checks, unit-test record, automated consistency checks, and final audit |
| Figure outputs | `figures/` | Audited raster and vector figures |

## Reused public sources

1. Zhang et al. (2022), supplementary workbook `41467_2022_33606_MOESM4_ESM.xlsx`, Fig. 2C. The source file is associated with Nature Communications article DOI `10.1038/s41467-022-33606-0`. The repository records the local source checksum but does not redistribute the workbook.
2. Zhou et al. (2023), GEO accession `GSE168182`, file `GSE168182_yZSJ025.fa`. The repository records the source checksum but does not redistribute the GEO FASTA.
3. Saccharomyces Genome Database feature and phenotype records, represented here by the processed files in `data/processed/` and the versioned Essentiality Gate membership table.
4. Synthetic-chromosome GenBank accessions `CP013608`, `KC880027`, `CP135953`, and `JN020955`, included under `outputs/revised_model/input_sources/` for coordinate mapping.

The exact checksums and asserted record counts are in `outputs/revised_model/inputs/input_provenance.json`.

## Deliberate exclusions

- SRA, GEO, and FASTQ raw sequencing files.
- The third-party Zhang supplementary workbook and Zhou FASTA.
- Original per-scenario event logs, full formal-run trajectory tables, and per-scenario SQLite endpoint catalogs. The compact derived cascade-replay tables are included because they are the direct source data for the new analysis.
- Invalidated legacy model outputs and earlier fixed-coordinate analyses.
- Full manuscript Word files and local user paths; the standalone Section 3.8 figure/table placement package is included as an analysis artifact.

## Reproducibility boundary

The included `endpoint_counts.csv.gz` is sufficient to reproduce the main diversity calculations, and `predicted_lox_frequency.csv` is sufficient to reproduce lox-frequency backtesting and null controls. The included cascade summary and compressed replay tables support verification and figure regeneration. Repeating the cascade replay from the original accepted-event stream requires running the formal workflow to regenerate the excluded scenario-level event logs. Recomputing trajectory-level gene-risk counts or endpoint-isomorphism mapping from raw simulated trajectories likewise requires those regenerated artifacts.

## Archival status

This GitHub repository is a version-controlled code and source-data location. A persistent dataset DOI has not been assigned in this repository. Before journal submission, archive a tagged release in an appropriate repository such as Zenodo and cite the resulting DOI. Author confirmation is required for the final data licence and creator metadata.
