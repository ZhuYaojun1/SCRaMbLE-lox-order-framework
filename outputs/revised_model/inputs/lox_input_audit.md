# Revised lox input audit

The Zhang 2022 Fig. 2C parser now requires both a valid header order and a non-empty chromosome-specific value.
The 604 padded blank cells are excluded from the active topology, candidate-pair construction, and predicted-frequency output.

## Fig. 2C row audit

| chromosome | raw_header_columns | nonempty_records | removed_blank_records | unique_lox_id | duplicate_lox_id | lox_order_min | lox_order_max | candidate_same_chromosome_pairs | sites_with_observed_frequency | removed_blank_examples |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| SynII | 273 | 273 | 0 | 273 | 0 | 0 | 272 | 37128 | 273 |  |
| SynIII | 273 | 100 | 173 | 100 | 0 | 0 | 99 | 4950 | 100 | 100;101;102;103;104;105;106;107;108;109 |
| SynIXR | 273 | 44 | 229 | 44 | 0 | 0 | 43 | 946 | 44 | 44;45;46;47;48;49;50;51;52;53 |
| SynVI | 273 | 71 | 202 | 71 | 0 | 0 | 70 | 2485 | 71 | 71;72;73;74;75;76;77;78;79;80 |

## Physical-coordinate evidence

| chromosome | fasta_record | fasta_length | fig2c_nonempty_records | physical_loxpsym_hits | coordinate_complete | unresolved_records |
| --- | --- | --- | --- | --- | --- | --- |
| SynII | chr2 | 770032 | 273 | 273 | True | 0 |
| SynIII | chr3 | 272194 | 100 | 100 | True | 0 |
| SynVI | chr6 | 242745 | 71 | 71 | True | 0 |
| SynIXR | chr9 | 442490 | 44 | 43 | False | 1 |

The supplied local FASTA resolves 487 of 488 Fig. 2C records. SynIXR contains 44 non-empty experimental records but 43 verifiable loxPsym motifs; the unmatched record is retained as an order-defined topology node with unresolved bp coordinate and is excluded from gene-level risk mapping.
