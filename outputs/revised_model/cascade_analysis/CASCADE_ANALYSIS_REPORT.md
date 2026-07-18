# Rearrangement cascade analysis

## Principal finding

Among 438,708 Essentiality-Gate-failing trajectories, 196,246 (44.73%) contained at least one inversion or duplication before the terminal deletion. Exact state replay evaluated 328,098 such precursor events with no topology, copy-number, or Gate-status mismatches. Prior rearrangements therefore occurred in a substantial fraction of terminal deletion histories, although their effects were heterogeneous rather than uniformly risk enhancing.

The copy-number Gate makes deletion the only event capable of producing an immediate Gate failure. The analysis consequently asks a narrower and testable question: whether earlier inversion or duplication changed the accessibility, proposal probability, or Gate consequence of the deletion path that later ended the trajectory.

## Failure-history composition

| history_class | gate_failing_trajectories |
|---|---|
| deletion_only | 215447 |
| inversion_only_history | 94112 |
| earlier_deletion_with_nondeletion | 49483 |
| earlier_deletion_only | 27015 |
| mixed_inversion_duplication | 26561 |
| duplication_only_history | 26090 |

## Exact terminal-pair effects

| precursor_event_type | n_precursor_events | n_trajectories | final_pair_activation_fraction | proposal_probability_increase_fraction | final_pair_lethal_hazard_increase_fraction | median_proposal_log2_ratio_when_nonzero |
|---|---|---|---|---|---|---|
| duplication | 86890 | 74133 | 0.0396 | 0.0396 | 0.0376 | -0.0151 |
| inversion | 241208 | 162430 | 0.0000 | 0.0340 | 0.0366 | 0.0000 |

| precursor_event_type | n_trajectories | trajectory_fraction_with_any_final_pair_hazard_increase | trajectory_fraction_ci_low | trajectory_fraction_ci_high | trajectory_fraction_with_pair_activation |
|---|---|---|---|---|---|
| inversion | 162430 | 0.0529 | 0.0518 | 0.0540 | 0.0000 |
| duplication | 74133 | 0.0441 | 0.0427 | 0.0456 | 0.0464 |

Inversions most often left the exact future terminal pair unchanged, but a minority shortened its current lox-order distance, changed its sampling weight, or converted its interval into a Gate-failing deletion span. Duplications could create a new copy-specific terminal pair, but more often expanded the competing pair space or buffered essential-gene copy number.

## Large terminal-pair probability changes

| precursor_event_type | change_class | n_precursor_events | n_trajectories |
|---|---|---|---|
| inversion | at_least_2_fold | 1155 | 1153 |
| inversion | at_least_10_fold | 110 | 110 |
| inversion | at_least_100_fold | 7 | 7 |
| inversion | new_pair_activation | 0 | 0 |
| duplication | at_least_2_fold | 0 | 0 |
| duplication | at_least_10_fold | 0 | 0 |
| duplication | at_least_100_fold | 0 | 0 |
| duplication | new_pair_activation | 3437 | 3437 |

Finite fold-change classes exclude newly activated pairs whose pre-event probability was zero. New-pair activation is therefore reported separately.

## Prospective all-pair Gate-failing deletion hazard

| precursor_event_type | n_sampled_events | hazard_increase_events | hazard_decrease_events | hazard_unchanged_events | hazard_increase_fraction | hazard_decrease_fraction | median_hazard_log2_ratio_nonzero |
|---|---|---|---|---|---|---|---|
| duplication | 4385 | 146 | 4239 | 0 | 0.0333 | 0.9667 | -0.0139 |
| inversion | 4400 | 1833 | 707 | 1860 | 0.4166 | 0.1607 | 0.0000 |

These all-pair values come from a deterministic, scenario-stratified sample of up to 50 precursor events of each type per scenario. They enumerate every active deletion pair in the state immediately before and after the precursor event and therefore summarize the one-step Gate-failing deletion hazard, not only the pair eventually observed in that same trajectory.

## Mechanistic categories

| precursor_event_type | cascade_mechanism | n_events | n_trajectories | event_fraction_within_type |
|---|---|---|---|---|
| duplication | proposal_space_renormalization | 74352 | 64485 | 0.8557 |
| duplication | distance_lengthening | 8434 | 8000 | 0.0971 |
| duplication | final_pair_activation | 3437 | 3437 | 0.0396 |
| duplication | no_direct_effect_on_final_pair | 541 | 477 | 0.0062 |
| duplication | essential_copy_buffering | 126 | 125 | 0.0015 |
| inversion | no_direct_effect_on_final_pair | 224487 | 154441 | 0.9307 |
| inversion | distance_shortening | 8132 | 7910 | 0.0337 |
| inversion | distance_lengthening | 7879 | 7661 | 0.0327 |
| inversion | lethal_span_activation | 705 | 705 | 0.0029 |
| inversion | essential_copy_buffering | 5 | 5 | 0.0000 |

## Missing-gene associations

| sgd_gene_id | gene_name | precursor_event_type | n_trajectories | hazard_increase_events | hazard_increase_event_fraction |
|---|---|---|---|---|---|
| S000000274 | ALG14 | inversion | 27956 | 2329 | 0.0565 |
| S000000283 | RPG1 | inversion | 28000 | 2321 | 0.0562 |
| S000000284 | SEC18 | inversion | 28000 | 2321 | 0.0562 |
| S000000356 | SPP381 | inversion | 27637 | 2308 | 0.0566 |
| S000000357 | RIB7 | inversion | 27637 | 2308 | 0.0566 |
| S000000358 | RPB5 | inversion | 27637 | 2308 | 0.0566 |
| S000000359 | CNS1 | inversion | 27637 | 2308 | 0.0566 |
| S000000295 | TIM12 | inversion | 28423 | 2299 | 0.0549 |
| S000000259 | PRP6 | inversion | 27649 | 2294 | 0.0562 |
| S000000371 | POP7 | inversion | 27615 | 2284 | 0.0561 |
| S000000306 | EXO84 | inversion | 28493 | 2283 | 0.0542 |
| S000000291 | RFC5 | inversion | 28250 | 2282 | 0.0548 |
| S000000292 | POL30 | inversion | 28250 | 2282 | 0.0548 |
| S000000293 | YBR089W | inversion | 28250 | 2282 | 0.0548 |
| S000000253 | REB1 | inversion | 27719 | 2270 | 0.0556 |

Gene rows describe associations with the essential genes absent after the later terminal deletion. They do not show that the precursor event deleted those genes or identify a biological lethal mechanism outside the modeled Essentiality Gate.

## Interpretation boundary

The exact terminal-pair analysis is conditioned on trajectories that later failed the Gate and is therefore retrospective; it quantifies path modification but does not estimate a population causal effect. The sampled all-pair analysis provides a broader state-conditioned risk measure, but it remains a one-step model hazard under the fitted event and pair-sampling rules. The results support a rearrangement-cascade interpretation in which earlier events can reshape later deletion accessibility, while also showing that most individual precursor events do not increase the eventual terminal-pair hazard.
