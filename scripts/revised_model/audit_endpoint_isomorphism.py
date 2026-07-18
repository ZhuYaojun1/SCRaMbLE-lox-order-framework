from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sqlite3
import zlib
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


PART_PATTERN = re.compile(r"([SL])\(([^,]+),([^,]+),([+-]\d+)\)")
REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit copy-label inflation in structural endpoint signatures.")
    parser.add_argument("--output-root", type=Path, default=REPOSITORY_ROOT / "outputs" / "revised_model")
    parser.add_argument("--rarefaction-reps", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260718)
    return parser.parse_args()


def diversity(counts: Counter[str]) -> dict[str, float]:
    values = np.asarray(list(counts.values()), dtype=float)
    probabilities = values / values.sum()
    shannon = float(-(probabilities * np.log(probabilities)).sum())
    return {
        "unique_endpoints": len(counts),
        "shannon_entropy": shannon,
        "effective_shannon_diversity": math.exp(shannon),
        "inverse_simpson_diversity": float(1.0 / np.sum(probabilities**2)),
    }


def rarefy(counts: Counter[str], depth: int, reps: int, rng: np.random.Generator) -> dict[str, float]:
    values = np.asarray(list(counts.values()), dtype=np.int64)
    metrics = []
    for _ in range(reps):
        sampled = rng.multivariate_hypergeometric(values, depth)
        sampled = sampled[sampled > 0]
        probabilities = sampled / sampled.sum()
        shannon = float(-(probabilities * np.log(probabilities)).sum())
        metrics.append((len(sampled), math.exp(shannon), float(1.0 / np.sum(probabilities**2))))
    matrix = np.asarray(metrics, dtype=float)
    return {
        "rarefied_unique_mean": float(matrix[:, 0].mean()),
        "rarefied_unique_ci_low": float(np.quantile(matrix[:, 0], 0.025)),
        "rarefied_unique_ci_high": float(np.quantile(matrix[:, 0], 0.975)),
        "rarefied_effective_shannon_mean": float(matrix[:, 1].mean()),
        "rarefied_effective_shannon_ci_low": float(np.quantile(matrix[:, 1], 0.025)),
        "rarefied_effective_shannon_ci_high": float(np.quantile(matrix[:, 1], 0.975)),
        "rarefied_inverse_simpson_mean": float(matrix[:, 2].mean()),
    }


def topology_payload(payload: str) -> str:
    chromosome, remainder = payload.split("|PARTS:", 1)
    parts_text = remainder.split("|JUNCTIONS:", 1)[0]
    parts = PART_PATTERN.findall(parts_text)
    if not parts:
        raise ValueError(f"Could not parse structural payload for {chromosome}")
    # Copy IDs are deliberately omitted. Ordered template identity, orientation,
    # multiplicity, lox identity, and all adjacencies remain encoded by this list.
    normalized = ",".join(f"{kind}({template},{orientation})" for kind, template, _copy, orientation in parts)
    return f"{chromosome}|PARTS:{normalized}"


def topology_signature(canonical: str, references: dict[str, str]) -> tuple[str, str]:
    if not canonical.startswith("segment-topology-v1||"):
        raise ValueError("Unexpected structural endpoint schema")
    normalized_chromosomes = []
    for token in canonical.split("||")[1:]:
        chromosome = token.split("|", 1)[0]
        payload = references[chromosome] if "|REFERENCE:" in token else token
        normalized_chromosomes.append(topology_payload(payload))
    normalized = "unlabeled-linear-topology-v1||" + "||".join(sorted(normalized_chromosomes))
    return hashlib.sha256(normalized.encode("utf-8")).hexdigest(), normalized


def load_canonicals(database: Path, hashes: list[str]) -> dict[str, str]:
    connection = sqlite3.connect(database)
    try:
        result: dict[str, str] = {}
        for start in range(0, len(hashes), 800):
            batch = hashes[start : start + 800]
            placeholders = ",".join("?" for _ in batch)
            query = f"SELECT structural_hash, canonical_zlib FROM endpoint_structure WHERE structural_hash IN ({placeholders})"
            for endpoint_hash, blob in connection.execute(query, batch):
                result[endpoint_hash] = zlib.decompress(blob).decode("utf-8")
        return result
    finally:
        connection.close()


def main() -> None:
    args = parse_args()
    main_run = args.output_root / "main_run"
    analysis = args.output_root / "analysis"
    analysis.mkdir(parents=True, exist_ok=True)
    endpoint_counts = pd.read_csv(main_run / "endpoint_counts.csv.gz")
    endpoint_counts = endpoint_counts[endpoint_counts["endpoint_type"] == "structural"].copy()
    population = pd.read_csv(main_run / "population_summary.csv")
    common_depth = int(population["gate_passing_trajectories"].min())
    rng = np.random.default_rng(args.seed)
    rows: list[dict[str, Any]] = []
    mapping_rows: list[dict[str, Any]] = []

    for scenario_index in range(1, 89):
        scenario_dir = main_run / "scenarios" / f"scenario_{scenario_index:03d}"
        config = json.loads((scenario_dir / "run_config.json").read_text(encoding="utf-8"))
        scenario = config["scenarios"][0]
        model_key, alpha, p_event = scenario
        seed = int(pd.read_csv(scenario_dir / "population_summary.csv")["seed"].iloc[0])
        selected = endpoint_counts[
            (endpoint_counts["model_key"] == model_key)
            & (endpoint_counts["alpha"] == float(alpha))
            & (endpoint_counts["p_event"] == float(p_event))
            & (endpoint_counts["seed"] == seed)
        ]
        references = json.loads((scenario_dir / "reference_structure_catalog.json").read_text(encoding="utf-8"))["chromosomes"]
        canonicals = load_canonicals(scenario_dir / "endpoint_catalog.sqlite", selected["endpoint_hash"].astype(str).tolist())
        if len(canonicals) != len(selected):
            raise RuntimeError(f"Scenario {scenario_index}: missing canonical endpoint payloads")

        provenance_counts: Counter[str] = Counter()
        topology_counts: Counter[str] = Counter()
        topology_members: dict[str, list[str]] = {}
        for record in selected.to_dict("records"):
            endpoint_hash = str(record["endpoint_hash"])
            count = int(record["count"])
            topology_hash, _ = topology_signature(canonicals[endpoint_hash], references)
            provenance_counts[endpoint_hash] += count
            topology_counts[topology_hash] += count
            topology_members.setdefault(topology_hash, []).append(endpoint_hash)
            mapping_rows.append(
                {
                    "scenario_index": scenario_index,
                    "model_key": model_key,
                    "alpha": alpha,
                    "p_event": p_event,
                    "seed": seed,
                    "provenance_endpoint_hash": endpoint_hash,
                    "unlabeled_topology_hash": topology_hash,
                    "trajectory_count": count,
                }
            )

        current = diversity(provenance_counts)
        normalized = diversity(topology_counts)
        current_rare = rarefy(provenance_counts, common_depth, args.rarefaction_reps, rng)
        normalized_rare = rarefy(topology_counts, common_depth, args.rarefaction_reps, rng)
        merged_groups = [members for members in topology_members.values() if len(members) > 1]
        trajectories_in_merged_groups = sum(topology_counts[key] for key, members in topology_members.items() if len(members) > 1)
        rows.append(
            {
                "scenario_index": scenario_index,
                "model_key": model_key,
                "alpha": alpha,
                "p_event": p_event,
                "seed": seed,
                "gate_passing_trajectories": sum(provenance_counts.values()),
                "common_rarefied_depth": common_depth,
                "provenance_unique_endpoints": current["unique_endpoints"],
                "unlabeled_topology_unique_endpoints": normalized["unique_endpoints"],
                "unique_inflation_ratio": current["unique_endpoints"] / normalized["unique_endpoints"],
                "provenance_effective_shannon": current["effective_shannon_diversity"],
                "unlabeled_topology_effective_shannon": normalized["effective_shannon_diversity"],
                "effective_shannon_inflation_ratio": current["effective_shannon_diversity"] / normalized["effective_shannon_diversity"],
                "provenance_inverse_simpson": current["inverse_simpson_diversity"],
                "unlabeled_topology_inverse_simpson": normalized["inverse_simpson_diversity"],
                "provenance_rarefied_effective_shannon_mean": current_rare["rarefied_effective_shannon_mean"],
                "unlabeled_topology_rarefied_effective_shannon_mean": normalized_rare["rarefied_effective_shannon_mean"],
                "rarefied_effective_shannon_inflation_ratio": current_rare["rarefied_effective_shannon_mean"] / normalized_rare["rarefied_effective_shannon_mean"],
                "unlabeled_topology_rarefied_ci_low": normalized_rare["rarefied_effective_shannon_ci_low"],
                "unlabeled_topology_rarefied_ci_high": normalized_rare["rarefied_effective_shannon_ci_high"],
                "topology_groups_merging_multiple_provenance_ids": len(merged_groups),
                "max_provenance_ids_per_topology": max((len(members) for members in topology_members.values()), default=1),
                "fraction_trajectories_in_merged_topology_groups": trajectories_in_merged_groups / sum(provenance_counts.values()),
            }
        )

    audit = pd.DataFrame(rows)
    mapping = pd.DataFrame(mapping_rows)
    audit.to_csv(analysis / "endpoint_isomorphism_audit.csv", index=False)
    mapping.to_csv(analysis / "endpoint_provenance_to_unlabeled_topology.csv.gz", index=False, compression="gzip")

    worst = audit.sort_values("rarefied_effective_shannon_inflation_ratio", ascending=False).iloc[0]
    best_topology = audit.sort_values("unlabeled_topology_rarefied_effective_shannon_mean", ascending=False).iloc[0]
    lines = [
        "# Structural endpoint isomorphism and copy-identity audit",
        "",
        "## Finding",
        "",
        "The current `segment-topology-v1` endpoint is not graph-isomorphism normalized. It orders chromosomes, but retains event-derived segment and lox copy IDs in PARTS, JUNCTIONS, and ACTIVE_LOX. Consequently, it is a lineage/provenance-aware signature rather than an unlabeled structural endpoint signature.",
        "",
        "The audit derived an unlabeled linear-topology signature by removing copy-instance labels while retaining chromosome identity, ordered segment/lox template identities, orientation, multiplicity, and adjacency. No Monte Carlo trajectories were rerun.",
        "",
        "## Quantitative effect",
        "",
        f"Across 88 scenarios, the median full-sample unique-endpoint inflation ratio was {audit['unique_inflation_ratio'].median():.3f} (range {audit['unique_inflation_ratio'].min():.3f}-{audit['unique_inflation_ratio'].max():.3f}).",
        f"The median rarefied effective-Shannon inflation ratio was {audit['rarefied_effective_shannon_inflation_ratio'].median():.3f} (range {audit['rarefied_effective_shannon_inflation_ratio'].min():.3f}-{audit['rarefied_effective_shannon_inflation_ratio'].max():.3f}).",
        f"The largest rarefied inflation occurred for {worst['model_key']}, alpha={worst['alpha']}, p_event={worst['p_event']}: {worst['rarefied_effective_shannon_inflation_ratio']:.3f}-fold.",
        f"After copy-label removal, the highest rarefied effective Shannon diversity was {best_topology['unlabeled_topology_rarefied_effective_shannon_mean']:.2f} for {best_topology['model_key']}, alpha={best_topology['alpha']}, p_event={best_topology['p_event']}.",
        "",
        "## Interpretation",
        "",
        "Copy-instance identity is useful only for lineage provenance and replay. Because IDs encode simulation event number and offset rather than an experimentally observable barcode, they should not define the primary structural endpoint diversity. Manuscript structural-diversity claims should use the unlabeled topology metric; the provenance-aware metric may be retained as a supplemental lineage-history measure.",
        "",
        "This normalization is appropriate for the current linear chromosome representation. It does not attempt to identify biologically equivalent structures under unmodeled circular symmetry or sequence-level homology.",
    ]
    (analysis / "endpoint_isomorphism_audit.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "scenarios": len(audit),
                "median_unique_inflation": float(audit["unique_inflation_ratio"].median()),
                "median_rarefied_effective_shannon_inflation": float(audit["rarefied_effective_shannon_inflation_ratio"].median()),
                "maximum_rarefied_effective_shannon_inflation": float(audit["rarefied_effective_shannon_inflation_ratio"].max()),
                "best_unlabeled_topology": best_topology.to_dict(),
            },
            indent=2,
            default=str,
        )
    )


if __name__ == "__main__":
    main()
