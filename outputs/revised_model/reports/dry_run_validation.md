# Dry-run validation

Overall status: **PASS**

| Check | Status | Value | Expected | Notes |
| --- | --- | --- | --- | --- |
| validated_lox_count | PASS | 488 | 488 |  |
| chromosome_lox_counts | PASS | {'SynII': 273, 'SynIII': 100, 'SynIXR': 44, 'SynVI': 71} | {'SynII': 273, 'SynIII': 100, 'SynIXR': 44, 'SynVI': 71} |  |
| candidate_pair_count | PASS | 45509 | 45509 |  |
| coordinate_resolved_lox | PASS | 487 | 487 |  |
| removed_padded_blank_records | PASS | 604 | 604 |  |
| predicted_rows_per_scenario | PASS | [488] | [488] |  |
| non_synii_direct_hic_pairs | PASS | 0 | 0 |  |
| reference_active_lox | PASS | 488 | 488 |  |
| reference_segment_count | PASS | 492 | 492 |  |
| accepted_event_reconciliation | PASS | 1531 | 1531 |  |
| trajectory_record_reconciliation | PASS | 800 | 800 |  |
| terminal_classification_reconciliation | PASS | 800 | 800 |  |
| endpoint_count_reconciliation | PASS | all endpoint types checked | each equals gate-passing trajectories |  |
| endpoint_catalog_reconstructable | PASS | sampled=20 | all canonical strings decompress |  |
| event_unit_tests | PASS | ./outputs\revised_model\unit_test_results.txt | seven tests pass |  |

## Coordinate scope

SynII, SynIII and SynVI have complete motif-supported lox coordinates for their Fig. 2C records. SynIXR has 44 non-empty Fig. 2C records but 43 verifiable local-FASTA loxPsym motifs. The unmatched record is retained only as an order-defined topology node; SynIXR is excluded from gene-level risk mapping.

## Event semantics

Deletion removes the current interval and one boundary copy; inversion reverses current segment order and orientation; duplication creates new segment and lox copy instances. Pair eligibility is rebuilt from the current topology at every accepted event.

Formal execution is permitted only when every automated check above passes.
