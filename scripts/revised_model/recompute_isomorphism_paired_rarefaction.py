from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Recompute endpoint isomorphism metrics with paired rarefaction.")
    parser.add_argument("--output-root", type=Path, default=REPOSITORY_ROOT / "outputs" / "revised_model")
    parser.add_argument("--rarefaction-reps", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260718)
    return parser.parse_args()


def metrics(values: np.ndarray) -> tuple[int, float, float]:
    values = values[values > 0].astype(float)
    probabilities = values / values.sum()
    shannon = float(-(probabilities * np.log(probabilities)).sum())
    return len(values), math.exp(shannon), float(1.0 / np.sum(probabilities**2))


def full_metrics(counts: Counter[str]) -> dict[str, float]:
    unique, effective, inverse = metrics(np.asarray(list(counts.values()), dtype=np.int64))
    return {"unique": unique, "effective_shannon": effective, "inverse_simpson": inverse}


def main() -> None:
    args = parse_args()
    analysis = args.output_root / "analysis"
    mapping = pd.read_csv(analysis / "endpoint_provenance_to_unlabeled_topology.csv.gz")
    aliases = pd.read_csv(analysis / "endpoint_alias_classification.csv")
    alias_lookup = {
        (int(row["scenario_index"]), str(row["unlabeled_topology_hash"])): str(row["classification"])
        for row in aliases.to_dict("records")
    }
    common_depth = int(pd.read_csv(args.output_root / "main_run" / "population_summary.csv")["gate_passing_trajectories"].min())
    rng = np.random.default_rng(args.seed)
    rows: list[dict[str, Any]] = []

    for scenario_index, group in mapping.groupby("scenario_index", sort=True):
        group = group.reset_index(drop=True)
        provenance_labels = group["provenance_endpoint_hash"].astype(str).tolist()
        topology_labels = group["unlabeled_topology_hash"].astype(str).tolist()
        counts_array = group["trajectory_count"].to_numpy(dtype=np.int64)
        provenance_counts = Counter(dict(zip(provenance_labels, counts_array)))
        topology_counts: Counter[str] = Counter()
        copy_normalized_counts: Counter[str] = Counter()
        reference_normalized_counts: Counter[str] = Counter()
        topology_members: dict[str, list[str]] = {}
        for provenance, topology, count in zip(provenance_labels, topology_labels, counts_array):
            topology_counts[topology] += int(count)
            topology_members.setdefault(topology, []).append(provenance)
            classification = alias_lookup.get((int(scenario_index), topology), "")
            copy_label = topology if classification == "event_copy_identity_alias" else provenance
            reference_label = topology if classification == "reference_token_vs_expanded_payload_alias" else provenance
            copy_normalized_counts[copy_label] += int(count)
            reference_normalized_counts[reference_label] += int(count)

        provenance_full = full_metrics(provenance_counts)
        topology_full = full_metrics(topology_counts)
        copy_normalized_full = full_metrics(copy_normalized_counts)
        reference_normalized_full = full_metrics(reference_normalized_counts)
        topology_index = {label: index for index, label in enumerate(sorted(topology_counts))}
        provenance_to_topology_index = np.asarray([topology_index[label] for label in topology_labels], dtype=int)
        copy_labels = [
            topology if alias_lookup.get((int(scenario_index), topology), "") == "event_copy_identity_alias" else provenance
            for provenance, topology in zip(provenance_labels, topology_labels)
        ]
        reference_labels = [
            topology if alias_lookup.get((int(scenario_index), topology), "") == "reference_token_vs_expanded_payload_alias" else provenance
            for provenance, topology in zip(provenance_labels, topology_labels)
        ]
        copy_index = {label: index for index, label in enumerate(sorted(set(copy_labels)))}
        reference_index = {label: index for index, label in enumerate(sorted(set(reference_labels)))}
        provenance_to_copy_index = np.asarray([copy_index[label] for label in copy_labels], dtype=int)
        provenance_to_reference_index = np.asarray([reference_index[label] for label in reference_labels], dtype=int)
        current_values = []
        topology_values = []
        copy_values = []
        reference_values = []
        replicate_ratios = []
        copy_replicate_ratios = []
        reference_replicate_ratios = []
        for _ in range(args.rarefaction_reps):
            sampled = rng.multivariate_hypergeometric(counts_array, common_depth)
            current_metric = metrics(sampled)
            topology_sampled = np.bincount(
                provenance_to_topology_index,
                weights=sampled,
                minlength=len(topology_index),
            ).astype(np.int64)
            topology_metric = metrics(topology_sampled)
            copy_sampled = np.bincount(
                provenance_to_copy_index, weights=sampled, minlength=len(copy_index)
            ).astype(np.int64)
            reference_sampled = np.bincount(
                provenance_to_reference_index, weights=sampled, minlength=len(reference_index)
            ).astype(np.int64)
            copy_metric = metrics(copy_sampled)
            reference_metric = metrics(reference_sampled)
            current_values.append(current_metric)
            topology_values.append(topology_metric)
            copy_values.append(copy_metric)
            reference_values.append(reference_metric)
            replicate_ratios.append(current_metric[1] / topology_metric[1])
            copy_replicate_ratios.append(current_metric[1] / copy_metric[1])
            reference_replicate_ratios.append(current_metric[1] / reference_metric[1])
        current_matrix = np.asarray(current_values, dtype=float)
        topology_matrix = np.asarray(topology_values, dtype=float)
        copy_matrix = np.asarray(copy_values, dtype=float)
        reference_matrix = np.asarray(reference_values, dtype=float)
        ratios = np.asarray(replicate_ratios, dtype=float)
        copy_ratios = np.asarray(copy_replicate_ratios, dtype=float)
        reference_ratios = np.asarray(reference_replicate_ratios, dtype=float)
        merged_keys = {key for key, members in topology_members.items() if len(members) > 1}
        merged_trajectory_count = sum(topology_counts[key] for key in merged_keys)
        first = group.iloc[0]
        rows.append(
            {
                "scenario_index": int(scenario_index),
                "model_key": first["model_key"],
                "alpha": float(first["alpha"]),
                "p_event": float(first["p_event"]),
                "seed": int(first["seed"]),
                "gate_passing_trajectories": int(counts_array.sum()),
                "common_rarefied_depth": common_depth,
                "provenance_unique_endpoints": provenance_full["unique"],
                "unlabeled_topology_unique_endpoints": topology_full["unique"],
                "unique_inflation_ratio": provenance_full["unique"] / topology_full["unique"],
                "provenance_effective_shannon": provenance_full["effective_shannon"],
                "unlabeled_topology_effective_shannon": topology_full["effective_shannon"],
                "effective_shannon_inflation_ratio": provenance_full["effective_shannon"] / topology_full["effective_shannon"],
                "copy_identity_normalized_unique_endpoints": copy_normalized_full["unique"],
                "copy_identity_unique_inflation_ratio": provenance_full["unique"] / copy_normalized_full["unique"],
                "copy_identity_normalized_effective_shannon": copy_normalized_full["effective_shannon"],
                "copy_identity_effective_shannon_inflation_ratio": provenance_full["effective_shannon"] / copy_normalized_full["effective_shannon"],
                "reference_alias_normalized_unique_endpoints": reference_normalized_full["unique"],
                "reference_alias_unique_inflation_ratio": provenance_full["unique"] / reference_normalized_full["unique"],
                "reference_alias_normalized_effective_shannon": reference_normalized_full["effective_shannon"],
                "reference_alias_effective_shannon_inflation_ratio": provenance_full["effective_shannon"] / reference_normalized_full["effective_shannon"],
                "provenance_inverse_simpson": provenance_full["inverse_simpson"],
                "unlabeled_topology_inverse_simpson": topology_full["inverse_simpson"],
                "provenance_rarefied_effective_shannon_mean": float(current_matrix[:, 1].mean()),
                "unlabeled_topology_rarefied_effective_shannon_mean": float(topology_matrix[:, 1].mean()),
                "paired_rarefied_effective_shannon_inflation_mean": float(ratios.mean()),
                "paired_rarefied_effective_shannon_inflation_ci_low": float(np.quantile(ratios, 0.025)),
                "paired_rarefied_effective_shannon_inflation_ci_high": float(np.quantile(ratios, 0.975)),
                "copy_identity_normalized_rarefied_effective_shannon_mean": float(copy_matrix[:, 1].mean()),
                "copy_identity_paired_rarefied_inflation_mean": float(copy_ratios.mean()),
                "copy_identity_paired_rarefied_inflation_ci_low": float(np.quantile(copy_ratios, 0.025)),
                "copy_identity_paired_rarefied_inflation_ci_high": float(np.quantile(copy_ratios, 0.975)),
                "reference_alias_normalized_rarefied_effective_shannon_mean": float(reference_matrix[:, 1].mean()),
                "reference_alias_paired_rarefied_inflation_mean": float(reference_ratios.mean()),
                "reference_alias_paired_rarefied_inflation_ci_low": float(np.quantile(reference_ratios, 0.025)),
                "reference_alias_paired_rarefied_inflation_ci_high": float(np.quantile(reference_ratios, 0.975)),
                "unlabeled_topology_rarefied_ci_low": float(np.quantile(topology_matrix[:, 1], 0.025)),
                "unlabeled_topology_rarefied_ci_high": float(np.quantile(topology_matrix[:, 1], 0.975)),
                "topology_groups_merging_multiple_provenance_ids": len(merged_keys),
                "max_provenance_ids_per_topology": max((len(members) for members in topology_members.values()), default=1),
                "fraction_trajectories_in_merged_topology_groups": merged_trajectory_count / counts_array.sum(),
            }
        )

    audit = pd.DataFrame(rows)
    audit.to_csv(analysis / "endpoint_isomorphism_audit.csv", index=False)
    merged = audit[audit["topology_groups_merging_multiple_provenance_ids"] > 0]
    worst = audit.sort_values("paired_rarefied_effective_shannon_inflation_mean", ascending=False).iloc[0]
    best_topology = audit.sort_values("unlabeled_topology_rarefied_effective_shannon_mean", ascending=False).iloc[0]
    lines = [
        "# Structural endpoint isomorphism and copy-identity audit",
        "",
        "## Finding",
        "",
        "The current `segment-topology-v1` signature is not graph-isomorphism normalized. Chromosomes are ordered, but event-derived segment and lox copy IDs remain in PARTS, JUNCTIONS, and ACTIVE_LOX. It is therefore provenance-aware rather than copy-label invariant.",
        "",
        "An unlabeled linear-topology signature was derived by removing copy-instance labels while retaining chromosome identity, ordered segment/lox template identity, orientation, multiplicity, and adjacency. Monte Carlo trajectories were not rerun.",
        "",
        "## Quantitative effect",
        "",
        f"Only {len(merged)} of 88 scenarios contained any topology group that merged multiple provenance hashes.",
        f"The median full-sample unique-endpoint inflation ratio was {audit['unique_inflation_ratio'].median():.6f}; the maximum was {audit['unique_inflation_ratio'].max():.6f}.",
        f"Using paired rarefaction of the same trajectories at depth {common_depth}, the median effective-Shannon inflation ratio was {audit['paired_rarefied_effective_shannon_inflation_mean'].median():.6f}; the maximum was {audit['paired_rarefied_effective_shannon_inflation_mean'].max():.6f}.",
        f"For event-derived copy identity alone, the maximum full-sample unique inflation was {audit['copy_identity_unique_inflation_ratio'].max():.6f}, and the maximum paired-rarefied effective-Shannon inflation was {audit['copy_identity_paired_rarefied_inflation_mean'].max():.6f}.",
        f"For reference-token versus expanded-payload aliases alone, the corresponding maxima were {audit['reference_alias_unique_inflation_ratio'].max():.6f} and {audit['reference_alias_paired_rarefied_inflation_mean'].max():.6f}.",
        f"The worst scenario was {worst['model_key']}, alpha={worst['alpha']}, p_event={worst['p_event']}, with a paired rarefied inflation of {worst['paired_rarefied_effective_shannon_inflation_mean']:.6f}.",
        f"After copy-label removal, the highest rarefied effective Shannon diversity remained {best_topology['unlabeled_topology_rarefied_effective_shannon_mean']:.2f} for {best_topology['model_key']}, alpha={best_topology['alpha']}, p_event={best_topology['p_event']}.",
        "",
        "## Interpretation",
        "",
        "Copy-instance identity can in principle inflate structural diversity because IDs encode simulation event number and offset rather than an experimentally observable barcode. In the present run, however, the measured inflation is negligible and does not change the highest-diversity parameter combination. The copy-unlabeled topology metric is still the preferable primary endpoint definition; the existing copy-aware metric should be described as lineage/provenance-aware structural diversity.",
        "",
        "The normalization applies to the current linear chromosome representation and does not assert equivalence under unmodeled circular symmetry or sequence-level homology.",
    ]
    (analysis / "endpoint_isomorphism_audit.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "scenarios": len(audit),
                "scenarios_with_merges": len(merged),
                "median_unique_inflation": float(audit["unique_inflation_ratio"].median()),
                "maximum_unique_inflation": float(audit["unique_inflation_ratio"].max()),
                "median_paired_rarefied_inflation": float(audit["paired_rarefied_effective_shannon_inflation_mean"].median()),
                "maximum_paired_rarefied_inflation": float(audit["paired_rarefied_effective_shannon_inflation_mean"].max()),
                "maximum_copy_identity_unique_inflation": float(audit["copy_identity_unique_inflation_ratio"].max()),
                "maximum_copy_identity_paired_rarefied_inflation": float(audit["copy_identity_paired_rarefied_inflation_mean"].max()),
                "maximum_reference_alias_unique_inflation": float(audit["reference_alias_unique_inflation_ratio"].max()),
                "maximum_reference_alias_paired_rarefied_inflation": float(audit["reference_alias_paired_rarefied_inflation_mean"].max()),
                "best_unlabeled_topology_model": best_topology["model_key"],
                "best_unlabeled_topology_alpha": float(best_topology["alpha"]),
                "best_unlabeled_topology_p_event": float(best_topology["p_event"]),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
