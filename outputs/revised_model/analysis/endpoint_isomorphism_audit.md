# Structural endpoint isomorphism and copy-identity audit

## Finding

The current `segment-topology-v1` signature is not graph-isomorphism normalized. Chromosomes are ordered, but event-derived segment and lox copy IDs remain in PARTS, JUNCTIONS, and ACTIVE_LOX. It is therefore provenance-aware rather than copy-label invariant.

An unlabeled linear-topology signature was derived by removing copy-instance labels while retaining chromosome identity, ordered segment/lox template identity, orientation, multiplicity, and adjacency. Monte Carlo trajectories were not rerun.

## Quantitative effect

Only 10 of 88 scenarios contained any topology group that merged multiple provenance hashes.
The median full-sample unique-endpoint inflation ratio was 1.000000; the maximum was 1.000536.
Using paired rarefaction of the same trajectories at depth 1144, the median effective-Shannon inflation ratio was 1.000000; the maximum was 1.001964.
For event-derived copy identity alone, the maximum full-sample unique inflation was 1.000428, and the maximum paired-rarefied effective-Shannon inflation was 1.001546.
For reference-token versus expanded-payload aliases alone, the corresponding maxima were 1.000357 and 1.001927.
The worst scenario was linear_distance, alpha=2.0, p_event=0.1, with a paired rarefied inflation of 1.001964.
After copy-label removal, the highest rarefied effective Shannon diversity remained 1138.78 for linear_distance, alpha=2.0, p_event=0.3.

## Interpretation

Copy-instance identity can in principle inflate structural diversity because IDs encode simulation event number and offset rather than an experimentally observable barcode. In the present run, however, the measured inflation is negligible and does not change the highest-diversity parameter combination. The copy-unlabeled topology metric is still the preferable primary endpoint definition; the existing copy-aware metric should be described as lineage/provenance-aware structural diversity.

The normalization applies to the current linear chromosome representation and does not assert equivalence under unmodeled circular symmetry or sequence-level homology.
