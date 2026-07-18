# Cascade analysis validation

Overall status: PASS

| Check | Status | Observed |
|---|---|---|
| All 88 scenarios were audited | PASS | 88 |
| Replay mismatches are zero | PASS | 0 |
| Precursor detail row count matches audit | PASS | 328098 |
| Precursor counts match audit scenario by scenario | PASS | matched |
| Trajectory detail row count matches audit | PASS | 196246 |
| Trajectory counts match audit scenario by scenario | PASS | matched |
| Failure-history categories close to all Gate failures | PASS | 438708 |
| Event summary and mechanism categories agree | PASS | {'duplication': 86890, 'inversion': 241208} |
| Total-hazard sample count matches audit | PASS | 8785 |
| Total-hazard direction categories are exhaustive | PASS | matched |
| All proposal probabilities are within [0, 1] | PASS | 0 |
| Pair-activation flags agree with pair existence | PASS | 0 |
| Lethality-activation flags agree with Gate status | PASS | 0 |

The validation checks accounting closure, scenario-level agreement, categorical exhaustiveness, probability bounds, and logical consistency of activation flags. It does not convert the retrospective terminal-pair analysis into a causal population estimate.
