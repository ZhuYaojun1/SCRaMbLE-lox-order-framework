# Data dictionary

## Principal audited inputs

| File | Unit of observation | Key fields |
|---|---|---|
| `lox_sites_validated.csv` | one active loxPsym site | chromosome, lox ID, lox order, observed count/frequency, physical coordinate, coordinate status |
| `lox_input_audit.csv` | one synthetic chromosome | raw columns, retained records, removed blanks, pair count, observed-frequency coverage |
| `hic_pair_coverage.csv` | one chromosome plus an all-chromosome total | candidate pairs, direct Hi-C pairs, fallback pairs, missing/zero/unresolved pairs |
| `hic_direct_pair_weights.csv` | one directly matched SynII lox pair | lox IDs, Hi-C bins, direct contact weight |
| `gene_lox_exact_mapping.csv` | one mapped ORF | SGD gene ID, synthetic coordinate, flanking lox IDs, segment ID, mapping status, risk eligibility |
| `reference_segments.csv` | one reference segment | chromosome, segment ID, flanking lox identities, ORFs, Essentiality Gate ORFs, structural features |
| `survival_gate_gene_sets.csv` | one gene-by-gate membership | gate name, SGD gene ID, chromosome, evidence, selection rule, gate version |

`position_bp`, feature coordinates, and segment lengths are measured in base pairs. `lox_order` and topological pair separation are ordinal values, not base-pair distances.

## Main-grid source data

| File | Unit of observation | Key fields |
|---|---|---|
| `population_summary.csv` | one parameter scenario | model, alpha, p_event, seed, initialized/passing/failing trajectories, accepted events, endpoint summaries |
| `endpoint_counts.csv.gz` | one endpoint state within a scenario | scenario, endpoint type, endpoint hash, trajectory count |
| `predicted_lox_frequency.csv` | one lox site within a scenario | scenario, model, alpha, p_event, chromosome, lox ID, predicted and observed normalized frequency |
| `scenario_execution_manifest.csv` | one formal scenario | return code, elapsed time, relative log location, reproducible command |

The main grid contains 88 scenarios: two primary sampling models, four alpha values, and eleven p_event values. The uniform-random implementation is an auxiliary baseline and is not part of the 88-scenario main grid.

## Analysis outputs

| File | Interpretation |
|---|---|
| `rarefied_diversity.csv` | 100 repeated without-replacement rarefactions at common depth 1,144, reported separately for structural and ORF copy-number endpoints |
| `full_sample_diversity.csv` | unstandardized endpoint diversity at each scenario's full passing-trajectory depth |
| `diversity_comparison.csv` | paired full-sample and rarefied metrics |
| `survival_diversity_pareto_frontier.csv` | scenarios evaluated on gate-passing fraction and rarefied structural diversity |
| `backtesting_report.csv` | Pearson, Spearman, RMSE, hotspot AUROC/AUPRC, and aligned-site counts |
| `null_control_summary.csv` | observed statistic and empirical one-sided probability for four 10,000-replicate controls |
| `chromosome_wise_holdout.csv` | chromosome-held-out parameter-selection performance |
| `exact_coordinate_gene_risk.csv` | strict-gate genes observed in gate-failing trajectories, limited to exact-mapped chromosomes |
| `endpoint_isomorphism_audit.csv` | provenance-aware versus copy-label-invariant endpoint comparisons |

## Terminology

- **Essentiality Gate-passing fraction:** fraction of initialized trajectories retaining at least one copy of every gene in the selected gate.
- **Gate-failing trajectory:** a model state with zero copies of at least one selected gate gene; this is not an experimentally observed dead-cell genome.
- **Structural endpoint:** canonical ordered segment/lox topology with orientation, adjacency, and active-copy identity.
- **ORF copy-number endpoint:** gene-copy vector reported separately from structural topology.
- **Partially Hi-C-informed sampling with distance-based fallback:** direct SynII contact mass plus current-topology distance fallback; non-SynII and unmatched pairs use fallback only.
