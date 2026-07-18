from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import sqlite3
import time
import zlib
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, TextIO

import numpy as np
import pandas as pd

from common import sha256_file
from segment_engine import PairSampler, clone_reference, count_parts, load_reference_genome


DEFAULT_P_EVENTS = [0.03, 0.05, 0.07, 0.09, 0.10, 0.12, 0.15, 0.18, 0.20, 0.25, 0.30]
DEFAULT_ALPHAS = [0.5, 1.0, 1.5, 2.0]
DISPLAY_NAMES = {
    "linear_distance": "linear-distance sampling",
    "partial_hic_fallback": "partially Hi-C-informed sampling with distance-based fallback",
    "uniform_random": "uniform-random auxiliary baseline",
}
REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the revised segment-level SCRaMbLE simulation.")
    parser.add_argument("--project-root", type=Path, default=REPOSITORY_ROOT)
    parser.add_argument("--output-root", type=Path, default=REPOSITORY_ROOT / "outputs" / "revised_model")
    parser.add_argument("--n-trajectories", type=int, default=10_000)
    parser.add_argument("--n-steps", type=int, default=20)
    parser.add_argument("--seed", type=int, default=202606)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--include-uniform", action="store_true")
    parser.add_argument("--scenario-index", type=int, default=None, help="Run one 1-based scenario from the formal grid.")
    parser.add_argument("--reuse-gate-files", action="store_true")
    parser.add_argument("--model-key", choices=sorted(DISPLAY_NAMES), default=None)
    parser.add_argument("--alpha", type=float, default=None)
    parser.add_argument("--p-event", type=float, default=None)
    parser.add_argument("--gate-name", default="strict")
    parser.add_argument("--run-subdir", default=None)
    return parser.parse_args()


def build_gate_sets(output_root: Path, project_root: Path) -> dict[str, set[str]]:
    inputs = output_root / "inputs"
    genes = pd.read_csv(inputs / "gene_lox_exact_mapping.csv", low_memory=False)
    candidates = pd.read_csv(
        project_root / "data" / "processed" / "essential_gene_candidates.csv", low_memory=False
    )
    eligible = genes[genes["eligible_for_gene_risk"].astype(str).str.lower().isin(["true", "1"])].copy()
    eligible = eligible[eligible["sgd_gene_id"].notna() & eligible["sgd_gene_id"].astype(str).ne("")]
    strict = set(eligible.loc[eligible["essential_status"] == "Essential", "sgd_gene_id"].astype(str))
    ambiguous = set(eligible.loc[eligible["essential_status"] == "Ambiguous", "sgd_gene_id"].astype(str))
    candidate_lookup = candidates.set_index("gene_id", drop=False).to_dict("index")
    high_confidence: set[str] = set()
    member_rows: list[dict[str, Any]] = []
    for gene_id in sorted(set(eligible["sgd_gene_id"].astype(str))):
        candidate = candidate_lookup.get(gene_id, {})
        phenotype = str(candidate.get("phenotype_summary", "")).lower()
        evidence_count = int(pd.to_numeric(pd.Series([candidate.get("evidence_count")]), errors="coerce").fillna(0).iloc[0])
        if evidence_count >= 1 and any(token in phenotype for token in ["inviable", "lethal"]):
            high_confidence.add(gene_id)
    gate_sets: dict[str, set[str]] = {
        "strict": strict,
        "high_confidence": high_confidence,
        "expanded": strict | ambiguous,
    }
    pool = sorted(set(eligible["sgd_gene_id"].astype(str)) - strict)
    rng = np.random.default_rng(20260717)
    for replicate in range(1, 6):
        size = min(len(strict), len(pool))
        gate_sets[f"random_control_{replicate}"] = set(rng.choice(pool, size=size, replace=False).tolist())

    eligible_lookup = eligible.set_index("sgd_gene_id", drop=False).to_dict("index")
    for gate_name, members in gate_sets.items():
        for gene_id in sorted(members):
            gene = eligible_lookup.get(gene_id, {})
            candidate = candidate_lookup.get(gene_id, {})
            if gate_name == "strict":
                rule = "exact-mapped ORF with phenotype-derived essential_status=Essential"
            elif gate_name == "expanded":
                rule = "strict set plus exact-mapped ORFs with essential_status=Ambiguous"
            elif gate_name == "high_confidence":
                rule = "exact-mapped ORF with phenotype summary containing inviable or lethal"
            else:
                rule = "size-matched random sample from exact-mapped non-strict ORFs"
            member_rows.append(
                {
                    "gate_name": gate_name,
                    "sgd_gene_id": gene_id,
                    "gene_name": gene.get("gene_name", ""),
                    "chromosome": gene.get("synthetic_chromosome", ""),
                    "essential_status": gene.get("essential_status", ""),
                    "evidence_count": candidate.get("evidence_count", ""),
                    "phenotype_summary": candidate.get("phenotype_summary", ""),
                    "selection_rule": rule,
                    "source_file": candidate.get("source_file", ""),
                    "gate_version": "revised-segment-model-v1",
                }
            )
    pd.DataFrame(member_rows).to_csv(inputs / "survival_gate_gene_sets.csv", index=False)
    summary = pd.DataFrame(
        [
            {
                "gate_name": name,
                "n_genes": len(members),
                "version": "revised-segment-model-v1",
                "conflict_rule": "Ambiguous excluded from strict and included only in expanded",
            }
            for name, members in gate_sets.items()
        ]
    )
    summary.to_csv(inputs / "survival_gate_summary.csv", index=False)
    return gate_sets


def load_gate_sets(output_root: Path) -> dict[str, set[str]]:
    members = pd.read_csv(output_root / "inputs" / "survival_gate_gene_sets.csv")
    return {
        name: set(group["sgd_gene_id"].astype(str))
        for name, group in members.groupby("gate_name")
    }


def event_probabilities(project_root: Path) -> dict[str, float]:
    frame = pd.read_csv(project_root / "data" / "processed" / "real_event_type_summary.csv")
    selected = frame[frame["condition"].astype(str).str.lower().eq("overall")]
    selected = selected[selected["event_type"].isin(["deletion", "inversion", "duplication"])]
    counts = selected.groupby("event_type")["event_count"].sum()
    total = float(counts.sum())
    if total <= 0:
        raise RuntimeError("No valid overall event-type counts.")
    return {event: float(counts[event] / total) for event in ["deletion", "inversion", "duplication"]}


class EndpointCatalog:
    def __init__(self, path: Path) -> None:
        self.connection = sqlite3.connect(path)
        self.connection.execute(
            "CREATE TABLE IF NOT EXISTS endpoint_structure "
            "(structural_hash TEXT PRIMARY KEY, canonical_zlib BLOB NOT NULL, schema_version TEXT NOT NULL)"
        )
        self.connection.execute(
            "CREATE TABLE IF NOT EXISTS endpoint_orf "
            "(orf_hash TEXT PRIMARY KEY, canonical_zlib BLOB NOT NULL, schema_version TEXT NOT NULL)"
        )

    def add(self, structural_hash: str, canonical: str, orf_hash: str, orf_canonical: str) -> None:
        self.connection.execute(
            "INSERT OR IGNORE INTO endpoint_structure VALUES (?, ?, ?)",
            (structural_hash, zlib.compress(canonical.encode("utf-8"), 9), "segment-topology-v1"),
        )
        self.connection.execute(
            "INSERT OR IGNORE INTO endpoint_orf VALUES (?, ?, ?)",
            (orf_hash, zlib.compress(orf_canonical.encode("utf-8"), 9), "orf-copy-number-v1"),
        )

    def commit(self) -> None:
        self.connection.commit()

    def close(self) -> None:
        self.commit()
        self.connection.close()


@dataclass
class ScenarioResult:
    summary: dict[str, Any]
    structural_counts: Counter[str]
    orf_counts: Counter[str]
    lox_counts_by_chromosome: dict[str, Counter[str]]


def diversity(counter: Counter[str]) -> dict[str, float]:
    total = sum(counter.values())
    if total == 0:
        return {"unique": 0, "shannon": 0.0, "effective_shannon": 0.0, "inverse_simpson": 0.0}
    probabilities = np.asarray(list(counter.values()), dtype=float) / total
    shannon = float(-(probabilities * np.log(probabilities)).sum())
    return {
        "unique": len(counter),
        "shannon": shannon,
        "effective_shannon": float(math.exp(shannon)),
        "inverse_simpson": float(1.0 / np.sum(probabilities**2)),
    }


def run_scenario(
    reference,
    sampler: PairSampler,
    gate_genes: set[str],
    event_probs: dict[str, float],
    n_trajectories: int,
    n_steps: int,
    p_event: float,
    seed: int,
    model_name: str,
    alpha: float,
    event_writer: csv.DictWriter,
    trajectory_writer: csv.DictWriter,
    catalog: EndpointCatalog,
    progress_prefix: str,
) -> ScenarioResult:
    rng = np.random.default_rng(seed)
    event_names = np.asarray(list(event_probs), dtype=object)
    event_weights = np.asarray(list(event_probs.values()), dtype=float)
    event_weights /= event_weights.sum()
    survivors = 0
    gate_failures = 0
    accepted_events = 0
    event_type_counts: Counter[str] = Counter()
    pair_source_counts: Counter[str] = Counter()
    structural_counts: Counter[str] = Counter()
    orf_counts: Counter[str] = Counter()
    lox_counts: dict[str, Counter[str]] = defaultdict(Counter)
    missing_counter: Counter[str] = Counter()
    start_time = time.perf_counter()

    for trajectory_id in range(n_trajectories):
        if trajectory_id and trajectory_id % max(100, n_trajectories // 10) == 0:
            print(
                f"{progress_prefix} {trajectory_id}/{n_trajectories}, elapsed={time.perf_counter() - start_time:.1f}s",
                flush=True,
            )
        state = clone_reference(reference)
        gate_status = "pass"
        terminal_step = n_steps
        missing: list[str] = []
        trajectory_event_count = 0
        for step in range(1, n_steps + 1):
            if rng.random() >= p_event:
                continue
            try:
                chromosome, left, right, pair_source = sampler.sample(state, rng)
            except RuntimeError:
                break
            event_type = str(rng.choice(event_names, p=event_weights))
            outcome = state.apply_event(event_type, chromosome, left, right)
            trajectory_event_count += 1
            accepted_events += 1
            event_type_counts[event_type] += 1
            pair_source_counts[pair_source] += 1
            lox_counts[chromosome][outcome.left_lox_template_id] += 1
            lox_counts[chromosome][outcome.right_lox_template_id] += 1
            missing = state.gate_missing(gate_genes)
            event_writer.writerow(
                {
                    "model_key": model_name,
                    "model_name": DISPLAY_NAMES[model_name],
                    "alpha": alpha,
                    "p_event": p_event,
                    "seed": seed,
                    "trajectory_id": trajectory_id,
                    "step": step,
                    "event_id": outcome.event_id,
                    "event_type": event_type,
                    "pair_weight_source": pair_source,
                    "chromosome": chromosome,
                    "left_lox_copy_id": outcome.left_lox_copy_id,
                    "right_lox_copy_id": outcome.right_lox_copy_id,
                    "left_lox_template_id": outcome.left_lox_template_id,
                    "right_lox_template_id": outcome.right_lox_template_id,
                    "selected_lox_separation": outcome.selected_lox_separation,
                    "affected_segment_copy_ids": ";".join(outcome.affected_segment_copy_ids),
                    "affected_orf_ids": ";".join(outcome.affected_orf_ids),
                    "affected_essential_orf_ids": ";".join(outcome.affected_essential_orf_ids),
                    "lox_count_before": outcome.lox_count_before,
                    "lox_count_after": outcome.lox_count_after,
                    "segment_count_before": outcome.segment_count_before,
                    "segment_count_after": outcome.segment_count_after,
                    "gate_status_after_event": "fail" if missing else "pass",
                    "missing_gate_genes": ";".join(missing),
                }
            )
            if missing:
                gate_status = "fail"
                terminal_step = step
                gate_failures += 1
                missing_counter.update(missing)
                break
        structural_hash, canonical = state.structural_signature()
        orf_hash, orf_canonical = state.orf_copy_number_signature()
        catalog.add(structural_hash, canonical, orf_hash, orf_canonical)
        if gate_status == "pass":
            survivors += 1
            structural_counts[structural_hash] += 1
            orf_counts[orf_hash] += 1
        segment_count, lox_count = count_parts(state)
        trajectory_writer.writerow(
            {
                "model_key": model_name,
                "model_name": DISPLAY_NAMES[model_name],
                "alpha": alpha,
                "p_event": p_event,
                "seed": seed,
                "trajectory_id": trajectory_id,
                "gate_status": gate_status,
                "terminal_step": terminal_step,
                "accepted_events": trajectory_event_count,
                "missing_gate_genes": ";".join(missing),
                "structural_endpoint_hash": structural_hash,
                "orf_copy_number_endpoint_hash": orf_hash,
                "final_segment_instances": segment_count,
                "final_active_lox_copies": lox_count,
            }
        )

    structural_diversity = diversity(structural_counts)
    orf_diversity = diversity(orf_counts)
    elapsed = time.perf_counter() - start_time
    summary = {
        "model_key": model_name,
        "model_name": DISPLAY_NAMES[model_name],
        "alpha": alpha,
        "p_event": p_event,
        "seed": seed,
        "initialized_trajectories": n_trajectories,
        "sequential_steps": n_steps,
        "accepted_events": accepted_events,
        "gate_passing_trajectories": survivors,
        "gate_failing_trajectories": gate_failures,
        "essentiality_gate_passing_fraction": survivors / n_trajectories,
        "essentiality_gate_failing_fraction": gate_failures / n_trajectories,
        "deletion_events": event_type_counts["deletion"],
        "inversion_events": event_type_counts["inversion"],
        "duplication_events": event_type_counts["duplication"],
        "direct_hic_events": pair_source_counts["direct_hic"],
        "distance_fallback_events": pair_source_counts["distance_fallback"],
        "uniform_events": pair_source_counts["uniform"],
        "structural_unique_endpoints_full_sample": structural_diversity["unique"],
        "structural_shannon_full_sample": structural_diversity["shannon"],
        "structural_effective_shannon_full_sample": structural_diversity["effective_shannon"],
        "structural_inverse_simpson_full_sample": structural_diversity["inverse_simpson"],
        "orf_unique_endpoints_full_sample": orf_diversity["unique"],
        "orf_shannon_full_sample": orf_diversity["shannon"],
        "orf_effective_shannon_full_sample": orf_diversity["effective_shannon"],
        "orf_inverse_simpson_full_sample": orf_diversity["inverse_simpson"],
        "top_missing_gate_genes": ";".join(f"{gene}:{count}" for gene, count in missing_counter.most_common(20)),
        "elapsed_seconds": elapsed,
    }
    return ScenarioResult(summary, structural_counts, orf_counts, dict(lox_counts))


def scenario_grid(dry_run: bool, include_uniform: bool) -> list[tuple[str, float, float]]:
    if dry_run:
        return [
            ("linear_distance", 1.0, 0.10),
            ("partial_hic_fallback", 1.0, 0.10),
            ("linear_distance", 2.0, 0.30),
            ("partial_hic_fallback", 0.5, 0.03),
        ]
    grid = [
        (model, alpha, p_event)
        for model in ["linear_distance", "partial_hic_fallback"]
        for alpha in DEFAULT_ALPHAS
        for p_event in DEFAULT_P_EVENTS
    ]
    if include_uniform:
        grid.extend(("uniform_random", 0.0, p_event) for p_event in DEFAULT_P_EVENTS)
    return grid


def main() -> None:
    args = parse_args()
    output = args.output_root
    inputs = output / "inputs"
    custom_scenario = args.model_key is not None or args.alpha is not None or args.p_event is not None
    if custom_scenario and not all(value is not None for value in [args.model_key, args.alpha, args.p_event]):
        raise ValueError("--model-key, --alpha, and --p-event must be supplied together.")
    if args.run_subdir:
        run_dir = output / args.run_subdir
    elif args.scenario_index is not None:
        run_dir = output / "main_run" / "scenarios" / f"scenario_{args.scenario_index:03d}"
    else:
        run_dir = output / ("dry_run" if args.dry_run else "main_run")
    run_dir.mkdir(parents=True, exist_ok=True)
    n_trajectories = min(args.n_trajectories, 200) if args.dry_run else args.n_trajectories
    reference = load_reference_genome(inputs)
    (run_dir / "reference_structure_catalog.json").write_text(
        json.dumps(
            {
                "schema": "lossless canonical chromosome payloads used by endpoint reference tokens",
                "chromosomes": reference.reference_payloads,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    (run_dir / "reference_orf_copy_catalog.json").write_text(
        json.dumps(
            {
                "schema": "reference ORF copy-number state",
                "orf_copy_numbers": dict(sorted(reference.orf_copy_numbers().items())),
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    gate_sets = load_gate_sets(output) if args.reuse_gate_files else build_gate_sets(output, args.project_root)
    if args.gate_name not in gate_sets:
        raise KeyError(f"Unknown gate set: {args.gate_name}")
    selected_gate = gate_sets[args.gate_name]
    reference_orfs = reference.orf_copy_numbers()
    missing_initial = sorted(gene for gene in selected_gate if reference_orfs.get(gene, 0) <= 0)
    if missing_initial:
        raise RuntimeError(f"Strict gate contains genes absent from the reference topology: {missing_initial[:20]}")
    event_probs = event_probabilities(args.project_root)
    direct = pd.read_csv(inputs / "hic_direct_pair_weights.csv")
    grid = scenario_grid(args.dry_run, args.include_uniform)
    if custom_scenario:
        grid = [(str(args.model_key), float(args.alpha), float(args.p_event))]
    elif args.scenario_index is not None:
        if args.dry_run:
            raise ValueError("--scenario-index is only valid for the formal grid.")
        if not 1 <= args.scenario_index <= len(grid):
            raise ValueError(f"Scenario index must be in 1..{len(grid)}")
        grid = [grid[args.scenario_index - 1]]

    event_fields = [
        "model_key", "model_name", "alpha", "p_event", "seed", "trajectory_id", "step", "event_id",
        "event_type", "pair_weight_source", "chromosome", "left_lox_copy_id", "right_lox_copy_id",
        "left_lox_template_id", "right_lox_template_id", "selected_lox_separation",
        "affected_segment_copy_ids", "affected_orf_ids", "affected_essential_orf_ids",
        "lox_count_before", "lox_count_after", "segment_count_before", "segment_count_after",
        "gate_status_after_event", "missing_gate_genes",
    ]
    trajectory_fields = [
        "model_key", "model_name", "alpha", "p_event", "seed", "trajectory_id", "gate_status",
        "terminal_step", "accepted_events", "missing_gate_genes", "structural_endpoint_hash",
        "orf_copy_number_endpoint_hash", "final_segment_instances", "final_active_lox_copies",
    ]
    event_handle = gzip.open(run_dir / "event_log.csv.gz", "wt", encoding="utf-8", newline="")
    trajectory_handle = gzip.open(run_dir / "trajectory_summary.csv.gz", "wt", encoding="utf-8", newline="")
    event_writer = csv.DictWriter(event_handle, fieldnames=event_fields)
    trajectory_writer = csv.DictWriter(trajectory_handle, fieldnames=trajectory_fields)
    event_writer.writeheader()
    trajectory_writer.writeheader()
    catalog = EndpointCatalog(run_dir / "endpoint_catalog.sqlite")
    summaries: list[dict[str, Any]] = []
    endpoint_rows: list[dict[str, Any]] = []
    frequency_rows: list[dict[str, Any]] = []
    lox_input = pd.read_csv(inputs / "lox_sites_validated.csv")
    lox_by_chrom = {
        chrom: group.sort_values("lox_order") for chrom, group in lox_input.groupby("chromosome")
    }

    try:
        for index, (model_name, alpha, p_event) in enumerate(grid, start=1):
            global_index = args.scenario_index if args.scenario_index is not None else index
            scenario_seed = args.seed if custom_scenario else int(
                np.random.SeedSequence([args.seed, global_index]).generate_state(1)[0]
            )
            sampler = PairSampler(model_name, alpha, direct if model_name == "partial_hic_fallback" else None)
            prefix = f"[scenario {global_index}] {model_name} alpha={alpha:g} p_event={p_event:g}"
            print(prefix, flush=True)
            result = run_scenario(
                reference, sampler, selected_gate, event_probs, n_trajectories, args.n_steps,
                p_event, scenario_seed, model_name, alpha, event_writer, trajectory_writer,
                catalog, prefix,
            )
            summaries.append(result.summary)
            for endpoint_type, counts in [("structural", result.structural_counts), ("orf_copy_number", result.orf_counts)]:
                for endpoint_hash, count in counts.items():
                    endpoint_rows.append(
                        {
                            "model_key": model_name,
                            "alpha": alpha,
                            "p_event": p_event,
                            "seed": scenario_seed,
                            "endpoint_type": endpoint_type,
                            "endpoint_hash": endpoint_hash,
                            "count": count,
                        }
                    )
            for chromosome, sites in lox_by_chrom.items():
                counts = result.lox_counts_by_chromosome.get(chromosome, Counter())
                chromosome_total = sum(counts.values())
                global_total = sum(sum(counter.values()) for counter in result.lox_counts_by_chromosome.values())
                for row in sites.to_dict("records"):
                    count = counts.get(str(row["lox_id"]), 0)
                    frequency_rows.append(
                        {
                            "model_key": model_name,
                            "model_name": DISPLAY_NAMES[model_name],
                            "alpha": alpha,
                            "p_event": p_event,
                            "seed": scenario_seed,
                            "chromosome": chromosome,
                            "lox_id": row["lox_id"],
                            "lox_order": row["lox_order"],
                            "position_bp": row.get("position_bp", ""),
                            "simulated_event_count": count,
                            "simulated_frequency_within_chromosome": count / chromosome_total if chromosome_total else 0.0,
                            "simulated_frequency_global": count / global_total if global_total else 0.0,
                            "observed_rearrangement_count": row["observed_rearrangement_count"],
                            "observed_frequency_within_chromosome": row["observed_rearrangement_frequency"],
                            "coordinate_status": row["coordinate_status"],
                        }
                    )
            catalog.commit()
            pd.DataFrame(summaries).to_csv(run_dir / "population_summary.partial.csv", index=False)
    finally:
        event_handle.close()
        trajectory_handle.close()
        catalog.close()

    pd.DataFrame(summaries).to_csv(run_dir / "population_summary.csv", index=False)
    pd.DataFrame(endpoint_rows).to_csv(run_dir / "endpoint_counts.csv.gz", index=False, compression="gzip")
    pd.DataFrame(frequency_rows).to_csv(run_dir / "predicted_lox_frequency.csv", index=False)
    config = {
        "run_type": "dry_run" if args.dry_run else "formal_main_grid",
        "base_seed": args.seed,
        "n_trajectories_per_scenario": n_trajectories,
        "n_steps": args.n_steps,
        "gate_name": args.gate_name,
        "gate_gene_count": len(selected_gate),
        "event_probabilities": event_probs,
        "scenario_count": len(grid),
        "global_scenario_index": args.scenario_index,
        "scenarios": grid,
        "terminology": {
            "survival_rate": "essentiality-gate passing fraction",
            "non_surviving": "gate-failing trajectory",
        },
        "input_hashes": {
            path.name: sha256_file(path)
            for path in [
                inputs / "lox_sites_validated.csv",
                inputs / "reference_segments.csv",
                inputs / "hic_direct_pair_weights.csv",
                inputs / "survival_gate_gene_sets.csv",
            ]
        },
    }
    (run_dir / "run_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")
    print(json.dumps({"run_dir": str(run_dir), "scenarios": len(grid), "summaries": len(summaries)}, indent=2))


if __name__ == "__main__":
    main()
