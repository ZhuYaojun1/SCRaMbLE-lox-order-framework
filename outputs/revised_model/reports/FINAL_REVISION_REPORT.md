# FINAL REVISION REPORT

## Scope and status

This report covers only the independent revised-model tree. Original raw data, legacy scripts, Stage 2.3 results, and manuscript files were not overwritten. The formal run uses 10,000 initialized trajectories per scenario and 20 sequential opportunities across 88 main-grid scenarios.

Automated checks: 17 passed and 0 failed.

## Original errors and corrections

1. **Padded lox columns.** The legacy Fig. 2C parser tested the shared order header but did not require a non-empty chromosome-specific count. It therefore emitted 604 blank padding cells as active sites. The revised parser requires a real value/record in the chromosome row. Active sites changed from 1,092 to 488 and same-chromosome candidate pairs from 148,512 to 45,509.
2. **Hi-C scope.** The legacy mapping made all 37,128 SynII pairs appear directly covered and could be described too broadly. Exact bp-to-bin mapping resolves 36,046 positive-contact SynII pairs; all non-SynII pairs use distance fallback. The model is named partially Hi-C-informed sampling with distance-based fallback.
3. **Gene-to-lox mapping.** Proportional midpoint scaling was removed. Features are mapped through accessioned synthetic-chromosome sequences and motif-supported lox coordinates. SynIXR has one unresolved physical lox record, so SynIXR gene-level risk ranking is disabled.
4. **State update engine.** Fixed-reference copy-number edits were replaced by ordered segment and lox copy instances. Deletion removes the interval and one recombination boundary; inversion reverses part order and orientation; duplication creates new segment and lox copy instances. Eligible pairs are rebuilt after every event.
5. **Endpoint definition.** The truncated changed-gene signature was removed. The provenance-aware endpoint losslessly encodes segment order, orientation, copy identity, junctions, and active lox identities. A separate copy-label-invariant linear-topology audit removes event-derived instance labels while retaining template order, orientation, multiplicity, and adjacency. ORF copy-number endpoints are reported separately.
6. **Diversity.** Main comparisons use 100 repeated without-replacement rarefactions to the minimum revised survivor depth. Full-sample results are supplemental. The discontinuous constrained-diversity score is retired in favor of a survival-diversity Pareto frontier and prespecified 0.50/0.70 retention indicators.
7. **Survival terminology.** Outputs use essentiality-gate passing/failing trajectory. The revised strict set contains only exact-mapped phenotype-derived Essential ORFs; it does not represent complete biological viability.

## Lox input audit

| chromosome | raw_header_columns | nonempty_records | removed_blank_records | unique_lox_id | lox_order_min | lox_order_max | candidate_same_chromosome_pairs | sites_with_observed_frequency |
|---|---|---|---|---|---|---|---|---|
| SynII | 273 | 273 | 0 | 273 | 0 | 272 | 37128 | 273 |
| SynIII | 273 | 100 | 173 | 100 | 0 | 99 | 4950 | 100 |
| SynIXR | 273 | 44 | 229 | 44 | 0 | 43 | 946 | 44 |
| SynVI | 273 | 71 | 202 | 71 | 0 | 70 | 2485 | 71 |

The 604 removed records do not appear in lox_sites_validated.csv, the dynamic active-lox topology, event proposal space, or predicted_lox_frequency.csv.

## Hi-C coverage audit

| chromosome | total_candidate_pairs | direct_hic_pairs | distance_fallback_pairs | unresolved_coordinate_pairs | direct_hic_fraction |
|---|---|---|---|---|---|
| SynII | 37128 | 36046 | 1082 | 0 | 0.9708575737987504 |
| SynIII | 4950 | 0 | 4950 | 0 | 0.0 |
| SynIXR | 946 | 0 | 946 | 43 | 0.0 |
| SynVI | 2485 | 0 | 2485 | 0 | 0.0 |
| ALL | 45509 | 36046 | 9463 | 43 | 0.7920631083961415 |

Direct and fallback proposal masses are normalized within the current eligible pair set. Direct contact weights are robust-scaled against the positive SynII contact distribution before mixture normalization; fallback remains the explicit topological-distance proposal rather than fabricated contact data.

## Exact coordinate mapping

| chromosome | source_accession | source_sequence_length | target_fasta_record | target_sequence_length | selected_features | mapped_features | unresolved_features |
|---|---|---|---|---|---|---|---|
| SynII | CP013608 | 769861 | chr2 | 770032 | 455 | 452 | 3 |
| SynIII | KC880027 | 272875 | chr3 | 272194 | 165 | 160 | 5 |
| SynVI | CP135953 | 249105 | chr6 | 242745 | 127 | 126 | 1 |
| SynIXR | JN020955 | 91010 | chr9 | 442490 | 55 | 54 | 1 |

Unresolved lox records: 1. Exact gene-level risk is limited to SynII, SynIII, and SynVI. A 20-ORF verification sample is saved in gene_lox_mapping_validation_sample.csv.

## Essentiality Gate

| gate_name | n_genes | version | conflict_rule |
|---|---|---|---|
| strict | 99 | revised-segment-model-v1 | Ambiguous excluded from strict and included only in expanded |
| high_confidence | 150 | revised-segment-model-v1 | Ambiguous excluded from strict and included only in expanded |
| expanded | 223 | revised-segment-model-v1 | Ambiguous excluded from strict and included only in expanded |
| random_control_1 | 99 | revised-segment-model-v1 | Ambiguous excluded from strict and included only in expanded |
| random_control_2 | 99 | revised-segment-model-v1 | Ambiguous excluded from strict and included only in expanded |
| random_control_3 | 99 | revised-segment-model-v1 | Ambiguous excluded from strict and included only in expanded |
| random_control_4 | 99 | revised-segment-model-v1 | Ambiguous excluded from strict and included only in expanded |
| random_control_5 | 99 | revised-segment-model-v1 | Ambiguous excluded from strict and included only in expanded |

The formal main grid uses the strict 99-gene set. Expanded, high-confidence, and five size-matched random controls are versioned independently. Gate failure means zero current copies of at least one selected ORF; it is not an experimentally observed dead-cell genome.

## Formal simulation results

The 88 scenarios initialized 880,000 trajectories and accepted 1,596,885 events: 706,106 deletions, 646,466 inversions, and 244,313 duplications.
The best rarefied structural effective Shannon diversity was 1138.15 for linear_distance at alpha=2.0 and p_event=0.3; its 95% rarefaction interval was 1130.46-1144.00.
At essentiality-gate passing fraction >=0.70, the largest rarefied structural diversity was 460.10 for linear_distance at alpha=2.0 and p_event=0.1.
The revised common rarefaction depth is 1144; all effective Shannon estimates satisfy D1 <= N. Full-versus-rarefied rank correlations are documented below.

| endpoint_type | spearman_full_vs_rarefied_rank | common_rarefied_depth |
|---|---|---|
| orf_copy_number | 0.993078792575112 | 1144 |
| structural | 0.8439470254658168 | 1144 |

## Endpoint isomorphism and copy-label audit

The original `segment-topology-v1` signature is not graph-isomorphism normalized because event-derived copy IDs remain encoded. After copy-label removal, 10 of 88 scenarios contained any merged provenance hashes. The maximum full-sample unique-endpoint inflation was 1.000536 and the maximum paired-rarefied effective-Shannon inflation was 1.001964. Event-derived copy identity alone produced maxima of 1.000428 for unique endpoints and 1.001546 for paired-rarefied effective Shannon diversity; reference-token representation aliases produced maxima of 1.000357 and 1.001927, respectively. The highest-diversity parameter combination remained linear-distance sampling, alpha=2.0, p_event=0.30. The copy-unlabeled topology metric is preferred for primary structural-diversity reporting; the copy-aware metric should be labeled lineage/provenance-aware.

## Frequency backtesting and controls

The largest partially Hi-C-informed Pearson correlation was r=0.4509 at alpha=2.0 and p_event=0.03. The best linear-distance Pearson correlation was r=0.4329 at alpha=1.0 and p_event=0.18.

| control_type | observed_pearson_r | replicates | null_mean | null_ci_low | null_ci_high | empirical_p_one_sided | model_key | alpha | p_event |
|---|---|---|---|---|---|---|---|---|---|
| global_site_label_shuffle | 0.4508658954193033 | 10000 | 0.0001456460885229 | -0.0799097449486242 | 0.1003755010795122 | 9.999000099990002e-05 | partial_hic_fallback | 2.0 | 0.03 |
| within_chromosome_site_label_permutation | 0.4508658954193033 | 10000 | 0.2652817810817974 | 0.1647882654479974 | 0.3756602108228716 | 0.0001999800019998 | partial_hic_fallback | 2.0 | 0.03 |
| between_chromosome_profile_permutation | 0.4508658954193033 | 10000 | 0.0642303536435243 | -0.185281002122393 | 0.492057264817794 | 0.081991800819918 | partial_hic_fallback | 2.0 | 0.03 |
| shuffled_hic_vs_distance_model_label | 0.059740689856653 | 10000 | 2.921030963637179e-05 | -0.1749131904930228 | 0.1787259322569638 | 0.2864713528647135 | partial_hic_fallback | 2.0 | 0.03 |

Chromosome-wise holdout contains 4 held-out chromosomes and remains a parameter-selection holdout within the same processed study, not independent experimental validation.

## Robustness analyses

The five-seed analysis contains 4 model-by-p_event summaries. Gate sensitivity contains 16 gate-by-p_event scenarios; complete member lists and rules are in survival_gate_gene_sets.csv.

## Exact-coordinate gene risk

The revised risk table contains 99 genes observed among gate-failing trajectories. All legacy named-gene rankings are invalid and must not be reused; only exact-coordinate records in exact_coordinate_gene_risk.csv may be discussed, with the gate-based limitation stated.

## Status of legacy results

### Invalidated

- The 1,092-site input, 148,512-pair proposal space, and 96,096 predicted-frequency row count.
- The claim that four chromosomes had direct Hi-C weighting and the old 37,128-direct-pair interpretation.
- All old effective/constrained diversity maxima, p_event optima, and threshold values derived from the fixed-reference engine.
- Pearson values 0.3022 and 0.319664 as evidence for the revised model; they belong to old analysis boundaries.
- The 133-ORF gate count, old PRP6/REB1/CHS2/CDS1/ALG14 risk ranking, and the old model-inferred dead-cell total.
- Any description of endpoint diversity as genome structural diversity when it was based only on truncated ORF copy-number states.

### Retained only at the conceptual level

- Survivor-only endpoint data motivate a risk-aware analysis of unobserved gate-failing trajectories.
- Increasing rearrangement opportunity creates a survival-diversity trade-off in the model, but its quantitative location must come from the revised tables.
- SynII Hi-C can be tested as a partial proposal constraint; any advantage must be reported from revised backtesting and null controls.

## Automated consistency checks

| check | status | observed | expected |
|---|---|---|---|
| active lox counts by chromosome | PASS | {"SynII": 273, "SynIII": 100, "SynIXR": 44, "SynVI": 71} | {"SynII": 273, "SynIII": 100, "SynVI": 71, "SynIXR": 44} |
| active lox total | PASS | 488 | 488 |
| removed blank padding | PASS | 604 | 604 |
| candidate pair total | PASS | 45509 | 45509 |
| direct Hi-C restricted to SynII | PASS | SynII | SynII only |
| formal scenario count | PASS | 88 | 88 |
| formal trajectory depth | PASS | 10000 | 10000 |
| trajectory balance | PASS | 88 | 88 |
| accepted event balance | PASS | 88 | 88 |
| predicted site rows | PASS | 42944 | 42944 |
| unique aligned sites per scenario | PASS | 488 | 488 |
| rarefied effective diversity mathematical bound | PASS | 176 | 176 |
| rarefaction replicates | PASS | 100 | >=100 |
| robustness return codes | PASS | 36 | 36 |
| endpoint isomorphism audit scenarios | PASS | 88 | 88 |
| unlabeled topology cannot increase unique endpoints | PASS | 88 | 88 |
| unlabeled rarefied diversity mathematical bound | PASS | 88 | 88 |

## Reproducibility boundary

The revised workflow does not claim a calibrated wet-lab induction protocol, complete fitness model, or base-pair breakpoint predictor. p_event remains a per-trajectory per-step model probability. All manuscript revision should be based only on this revised output tree and should preserve the distinction between essentiality-gate failure and biological death.
