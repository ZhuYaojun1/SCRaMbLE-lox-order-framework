# Public release checks

Validation date: 2026-07-19

| Check | Result | Observed |
|---|---|---|
| Python source compilation | PASS | all files under `scripts/revised_model/` compiled |
| Segment-engine unit tests | PASS | 7 of 7 tests |
| Public-path smoke simulation | PASS | 25 trajectories, 3 steps, 11 accepted events |
| Active lox records | PASS | 488 |
| Removed padded blank records | PASS | 604 |
| Candidate intrachromosomal pairs | PASS | 45,509 |
| Direct positive SynII Hi-C pairs | PASS | 36,046 |
| Main-grid scenarios | PASS | 88 |
| Predicted lox-frequency rows | PASS | 42,944 = 88 x 488 |
| Endpoint-count file readability | PASS | 478,374 aggregated rows |
| Cascade scenarios audited | PASS | 88 of 88 |
| Cascade replay accounting | PASS | 328,098 precursor events and 196,246 trajectories reconciled |
| Cascade replay mismatches | PASS | 0 |
| Cascade logical checks | PASS | 13 of 13 checks |
| Cascade manuscript source data rebuild | PASS | all six CSV source tables matched byte for byte |
| Cascade PNG rebuild | PASS | both 400-dpi PNG files matched byte for byte |
| Private absolute-path scan | PASS | no project-drive absolute path or named user path retained |
| GitHub single-file size limit | PASS | no file is 100 MB or larger |

The smoke simulation was executed in an ignored local validation directory and is not part of the deposited scientific results. The tracked main-grid and analysis files are copies of the audited revised-model outputs; only path strings in provenance fields were converted to repository-relative form.
