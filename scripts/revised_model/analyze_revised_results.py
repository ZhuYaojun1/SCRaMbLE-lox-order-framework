from __future__ import annotations

import argparse
import gzip
import itertools
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Analyze revised SCRaMbLE model outputs.")
    parser.add_argument("--output-root", type=Path, default=REPOSITORY_ROOT / "outputs" / "revised_model")
    parser.add_argument("--run-name", default="main_run")
    parser.add_argument("--rarefaction-reps", type=int, default=100)
    parser.add_argument("--null-reps", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=20260717)
    return parser.parse_args()


def correlation(x: np.ndarray, y: np.ndarray, rank: bool = False) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    x, y = x[valid], y[valid]
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    if rank:
        x = pd.Series(x).rank(method="average").to_numpy()
        y = pd.Series(y).rank(method="average").to_numpy()
    return float(np.corrcoef(x, y)[0, 1])


def auroc(labels: np.ndarray, scores: np.ndarray) -> float:
    labels = np.asarray(labels, dtype=int)
    scores = np.asarray(scores, dtype=float)
    positives = int(labels.sum())
    negatives = len(labels) - positives
    if positives == 0 or negatives == 0:
        return float("nan")
    ranks = pd.Series(scores).rank(method="average").to_numpy()
    return float((ranks[labels == 1].sum() - positives * (positives + 1) / 2) / (positives * negatives))


def average_precision(labels: np.ndarray, scores: np.ndarray) -> float:
    labels = np.asarray(labels, dtype=int)
    if labels.sum() == 0:
        return float("nan")
    order = np.argsort(-np.asarray(scores, dtype=float), kind="mergesort")
    sorted_labels = labels[order]
    precision = np.cumsum(sorted_labels) / np.arange(1, len(sorted_labels) + 1)
    return float(precision[sorted_labels == 1].mean())


def metric_record(group: pd.DataFrame) -> dict[str, float]:
    x = group["simulated_frequency_within_chromosome"].to_numpy(dtype=float)
    y = group["observed_frequency_within_chromosome"].to_numpy(dtype=float)
    threshold = group.groupby("chromosome")["observed_frequency_within_chromosome"].transform(
        lambda values: values.quantile(0.90)
    )
    labels = (group["observed_frequency_within_chromosome"] >= threshold).astype(int).to_numpy()
    return {
        "pearson_r": correlation(x, y),
        "spearman_r": correlation(x, y, rank=True),
        "rmse": float(np.sqrt(np.mean((x - y) ** 2))),
        "hotspot_auroc_top10_within_chromosome": auroc(labels, x),
        "hotspot_auprc_top10_within_chromosome": average_precision(labels, x),
        "hotspot_positive_fraction": float(labels.mean()),
        "n_sites": len(group),
    }


def full_sample_table(population: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for row in population.to_dict("records"):
        for endpoint_type, prefix in [("structural", "structural"), ("orf_copy_number", "orf")]:
            rows.append(
                {
                    "model_key": row["model_key"],
                    "model_name": row["model_name"],
                    "alpha": row["alpha"],
                    "p_event": row["p_event"],
                    "seed": row["seed"],
                    "endpoint_type": endpoint_type,
                    "sample_size": row["gate_passing_trajectories"],
                    "unique_endpoints": row[f"{prefix}_unique_endpoints_full_sample"],
                    "shannon_entropy": row[f"{prefix}_shannon_full_sample"],
                    "effective_shannon_diversity": row[f"{prefix}_effective_shannon_full_sample"],
                    "inverse_simpson_diversity": row[f"{prefix}_inverse_simpson_full_sample"],
                    "sample_size_standardized": False,
                }
            )
    return pd.DataFrame(rows)


def rarefy_counts(counts: np.ndarray, depth: int, rng: np.random.Generator) -> tuple[int, float, float, float]:
    sampled = rng.multivariate_hypergeometric(counts.astype(np.int64), depth)
    sampled = sampled[sampled > 0]
    probabilities = sampled / sampled.sum()
    shannon = float(-(probabilities * np.log(probabilities)).sum())
    inverse_simpson = float(1.0 / np.sum(probabilities**2))
    return len(sampled), shannon, float(math.exp(shannon)), inverse_simpson


def rarefaction(endpoint_counts: pd.DataFrame, population: pd.DataFrame, reps: int, seed: int) -> pd.DataFrame:
    min_survivors = int(population["gate_passing_trajectories"].min())
    if min_survivors <= 0:
        raise RuntimeError("Rarefaction requires at least one gate-passing endpoint in every scenario.")
    rng = np.random.default_rng(seed)
    rows: list[dict[str, Any]] = []
    keys = ["model_key", "alpha", "p_event", "seed", "endpoint_type"]
    for key, group in endpoint_counts.groupby(keys, sort=True):
        counts = group["count"].to_numpy(dtype=np.int64)
        total = int(counts.sum())
        if total < min_survivors:
            raise RuntimeError(f"Endpoint count below common depth for {key}: {total} < {min_survivors}")
        values = np.asarray([rarefy_counts(counts, min_survivors, rng) for _ in range(reps)], dtype=float)
        record: dict[str, Any] = dict(zip(keys, key))
        record.update(
            {
                "original_survivor_depth": total,
                "rarefied_survivor_depth": min_survivors,
                "rarefaction_replicates": reps,
            }
        )
        names = ["unique_endpoints", "shannon_entropy", "effective_shannon_diversity", "inverse_simpson_diversity"]
        for index, name in enumerate(names):
            record[f"{name}_mean"] = float(values[:, index].mean())
            record[f"{name}_ci_low"] = float(np.quantile(values[:, index], 0.025))
            record[f"{name}_ci_high"] = float(np.quantile(values[:, index], 0.975))
        rows.append(record)
    return pd.DataFrame(rows)


def diversity_comparison(full: pd.DataFrame, rare: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    merged = full.merge(
        rare,
        on=["model_key", "alpha", "p_event", "seed", "endpoint_type"],
        how="inner",
    )
    merged["full_sample_rank"] = merged.groupby("endpoint_type")["effective_shannon_diversity"].rank(
        ascending=False, method="min"
    )
    merged["rarefied_rank"] = merged.groupby("endpoint_type")["effective_shannon_diversity_mean"].rank(
        ascending=False, method="min"
    )
    rank_rows = []
    for endpoint_type, group in merged.groupby("endpoint_type"):
        rank_rows.append(
            {
                "endpoint_type": endpoint_type,
                "spearman_full_vs_rarefied_rank": correlation(
                    group["effective_shannon_diversity"].to_numpy(),
                    group["effective_shannon_diversity_mean"].to_numpy(),
                    rank=True,
                ),
                "common_rarefied_depth": int(group["rarefied_survivor_depth"].iloc[0]),
            }
        )
    return merged, pd.DataFrame(rank_rows)


def pareto_frontier(population: pd.DataFrame, rare: pd.DataFrame) -> pd.DataFrame:
    structural = rare[rare["endpoint_type"] == "structural"]
    merged = population.merge(
        structural[
            ["model_key", "alpha", "p_event", "seed", "effective_shannon_diversity_mean", "effective_shannon_diversity_ci_low", "effective_shannon_diversity_ci_high"]
        ],
        on=["model_key", "alpha", "p_event", "seed"],
    )
    passing = merged["essentiality_gate_passing_fraction"].to_numpy(dtype=float)
    diversity = merged["effective_shannon_diversity_mean"].to_numpy(dtype=float)
    pareto = []
    for index in range(len(merged)):
        dominated = np.any(
            (passing >= passing[index])
            & (diversity >= diversity[index])
            & ((passing > passing[index]) | (diversity > diversity[index]))
        )
        pareto.append(not bool(dominated))
    merged["pareto_optimal"] = pareto
    merged["passes_retention_0_50"] = merged["essentiality_gate_passing_fraction"] >= 0.50
    merged["passes_retention_0_70"] = merged["essentiality_gate_passing_fraction"] >= 0.70
    return merged.sort_values(["pareto_optimal", "essentiality_gate_passing_fraction"], ascending=[False, False])


def backtesting(predicted: pd.DataFrame) -> pd.DataFrame:
    keys = ["model_key", "model_name", "alpha", "p_event", "seed"]
    rows = []
    for key, group in predicted.groupby(keys, sort=True):
        record = dict(zip(keys, key))
        record.update(metric_record(group))
        rows.append(record)
    return pd.DataFrame(rows)


def holdout(predicted: pd.DataFrame) -> pd.DataFrame:
    scenario_keys = ["model_key", "model_name", "alpha", "p_event", "seed"]
    chromosomes = sorted(predicted["chromosome"].unique())
    rows: list[dict[str, Any]] = []
    for held_out in chromosomes:
        candidates: list[tuple[float, tuple[Any, ...], pd.DataFrame]] = []
        for key, group in predicted.groupby(scenario_keys, sort=True):
            training = group[group["chromosome"] != held_out]
            score = metric_record(training)["pearson_r"]
            candidates.append((score, key, group))
        candidates.sort(key=lambda item: (-np.inf if not np.isfinite(item[0]) else item[0]), reverse=True)
        training_score, key, selected = candidates[0]
        test = selected[selected["chromosome"] == held_out]
        record = dict(zip(scenario_keys, key))
        record.update(
            {
                "held_out_chromosome": held_out,
                "training_pearson_r": training_score,
                **{f"holdout_{name}": value for name, value in metric_record(test).items()},
            }
        )
        rows.append(record)
    return pd.DataFrame(rows)


def null_controls(predicted: pd.DataFrame, backtest: pd.DataFrame, reps: int, seed: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    hic = backtest[backtest["model_key"] == "partial_hic_fallback"].sort_values("pearson_r", ascending=False)
    selected_row = hic.iloc[0]
    selected = predicted[
        (predicted["model_key"] == selected_row["model_key"])
        & (predicted["alpha"] == selected_row["alpha"])
        & (predicted["p_event"] == selected_row["p_event"])
        & (predicted["seed"] == selected_row["seed"])
    ].copy()
    x = selected["simulated_frequency_within_chromosome"].to_numpy(dtype=float)
    y = selected["observed_frequency_within_chromosome"].to_numpy(dtype=float)
    observed = correlation(x, y)
    rng = np.random.default_rng(seed)
    global_values = np.empty(reps, dtype=float)
    within_values = np.empty(reps, dtype=float)
    chromosome_indices = [group.index.to_numpy() for _, group in selected.groupby("chromosome")]
    index_to_position = {index: position for position, index in enumerate(selected.index)}
    local_indices = [np.asarray([index_to_position[index] for index in group], dtype=int) for group in chromosome_indices]
    for replicate in range(reps):
        global_values[replicate] = correlation(x, rng.permutation(y))
        permuted = y.copy()
        for indices in local_indices:
            permuted[indices] = rng.permutation(permuted[indices])
        within_values[replicate] = correlation(x, permuted)
    rows = []
    distributions = []
    for control, values in [("global_site_label_shuffle", global_values), ("within_chromosome_site_label_permutation", within_values)]:
        empirical = (1 + int(np.sum(values >= observed))) / (len(values) + 1)
        rows.append(
            {
                "control_type": control,
                "observed_pearson_r": observed,
                "replicates": len(values),
                "null_mean": float(values.mean()),
                "null_ci_low": float(np.quantile(values, 0.025)),
                "null_ci_high": float(np.quantile(values, 0.975)),
                "empirical_p_one_sided": empirical,
                "model_key": selected_row["model_key"],
                "alpha": selected_row["alpha"],
                "p_event": selected_row["p_event"],
            }
        )
        distributions.extend(
            {"control_type": control, "replicate": index + 1, "pearson_r": value}
            for index, value in enumerate(values)
        )

    # Reassign complete within-chromosome profiles. Profiles are rank-preserving
    # interpolated because chromosome site counts differ.
    chromosomes = sorted(selected["chromosome"].unique())
    x_by_chromosome = {
        chromosome: selected.loc[selected["chromosome"] == chromosome, "simulated_frequency_within_chromosome"].to_numpy(dtype=float)
        for chromosome in chromosomes
    }
    y_by_chromosome = {
        chromosome: selected.loc[selected["chromosome"] == chromosome, "observed_frequency_within_chromosome"].to_numpy(dtype=float)
        for chromosome in chromosomes
    }

    def resize_profile(values: np.ndarray, size: int) -> np.ndarray:
        if len(values) == size:
            return values.copy()
        source_q = (np.arange(len(values)) + 0.5) / len(values)
        target_q = (np.arange(size) + 0.5) / size
        return np.interp(target_q, source_q, values)

    between_values = np.empty(reps, dtype=float)
    for replicate in range(reps):
        permuted_names = rng.permutation(chromosomes)
        x_blocks = []
        y_blocks = []
        for target, source in zip(chromosomes, permuted_names):
            x_blocks.append(x_by_chromosome[target])
            y_blocks.append(resize_profile(y_by_chromosome[str(source)], len(x_by_chromosome[target])))
        between_values[replicate] = correlation(np.concatenate(x_blocks), np.concatenate(y_blocks))
    rows.append(
        {
            "control_type": "between_chromosome_profile_permutation",
            "observed_pearson_r": observed,
            "replicates": reps,
            "null_mean": float(np.nanmean(between_values)),
            "null_ci_low": float(np.nanquantile(between_values, 0.025)),
            "null_ci_high": float(np.nanquantile(between_values, 0.975)),
            "empirical_p_one_sided": (1 + int(np.nansum(between_values >= observed))) / (reps + 1),
            "model_key": selected_row["model_key"],
            "alpha": selected_row["alpha"],
            "p_event": selected_row["p_event"],
        }
    )
    distributions.extend(
        {"control_type": "between_chromosome_profile_permutation", "replicate": index + 1, "pearson_r": value}
        for index, value in enumerate(between_values)
    )

    matched_linear = predicted[
        (predicted["model_key"] == "linear_distance")
        & (predicted["alpha"] == selected_row["alpha"])
        & (predicted["p_event"] == selected_row["p_event"])
    ].sort_values(["chromosome", "lox_order"])
    matched_hic = selected.sort_values(["chromosome", "lox_order"])
    if not matched_linear[["chromosome", "lox_id"]].reset_index(drop=True).equals(
        matched_hic[["chromosome", "lox_id"]].reset_index(drop=True)
    ):
        raise RuntimeError("Shuffled model-label control requires exact site alignment.")
    linear_x = matched_linear["simulated_frequency_within_chromosome"].to_numpy(dtype=float)
    hic_x = matched_hic["simulated_frequency_within_chromosome"].to_numpy(dtype=float)
    matched_y = matched_hic["observed_frequency_within_chromosome"].to_numpy(dtype=float)
    observed_advantage = correlation(hic_x, matched_y) - correlation(linear_x, matched_y)
    label_values = np.empty(reps, dtype=float)
    for replicate in range(reps):
        swap = rng.random(len(hic_x)) < 0.5
        permuted_hic = np.where(swap, linear_x, hic_x)
        permuted_linear = np.where(swap, hic_x, linear_x)
        label_values[replicate] = correlation(permuted_hic, matched_y) - correlation(permuted_linear, matched_y)
    rows.append(
        {
            "control_type": "shuffled_hic_vs_distance_model_label",
            "observed_pearson_r": observed_advantage,
            "replicates": reps,
            "null_mean": float(np.nanmean(label_values)),
            "null_ci_low": float(np.nanquantile(label_values, 0.025)),
            "null_ci_high": float(np.nanquantile(label_values, 0.975)),
            "empirical_p_one_sided": (1 + int(np.nansum(label_values >= observed_advantage))) / (reps + 1),
            "model_key": selected_row["model_key"],
            "alpha": selected_row["alpha"],
            "p_event": selected_row["p_event"],
        }
    )
    distributions.extend(
        {"control_type": "shuffled_hic_vs_distance_model_label", "replicate": index + 1, "pearson_r": value}
        for index, value in enumerate(label_values)
    )
    return pd.DataFrame(rows), pd.DataFrame(distributions)


def gene_risk(run_dir: Path, output_root: Path) -> pd.DataFrame:
    mapping = pd.read_csv(output_root / "inputs" / "gene_lox_exact_mapping.csv")
    mapping = mapping[mapping["eligible_for_gene_risk"].astype(str).str.lower().isin(["true", "1"])]
    lookup = mapping.set_index("sgd_gene_id", drop=False).to_dict("index")
    counts: Counter[str] = Counter()
    total_failures = 0
    trajectory_files = sorted((run_dir / "scenarios").glob("scenario_*/trajectory_summary.csv.gz"))
    if not trajectory_files and (run_dir / "trajectory_summary.csv.gz").exists():
        trajectory_files = [run_dir / "trajectory_summary.csv.gz"]
    for path in trajectory_files:
        for chunk in pd.read_csv(path, usecols=["gate_status", "missing_gate_genes"], chunksize=100_000):
            failed = chunk[chunk["gate_status"] == "fail"]
            total_failures += len(failed)
            for value in failed["missing_gate_genes"].dropna().astype(str):
                counts.update(token for token in value.split(";") if token)
    rows = []
    for gene_id, count in counts.most_common():
        record = lookup.get(gene_id, {})
        rows.append(
            {
                "sgd_gene_id": gene_id,
                "gene_name": record.get("gene_name", ""),
                "chromosome": record.get("synthetic_chromosome", ""),
                "left_lox_id": record.get("left_lox_id", ""),
                "right_lox_id": record.get("right_lox_id", ""),
                "segment_id": record.get("segment_id", ""),
                "gate_failure_count": count,
                "fraction_of_gate_failures": count / total_failures if total_failures else np.nan,
                "coordinate_scope": "exact-mapped chromosomes only; SynIXR excluded",
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    run_dir = args.output_root / args.run_name
    analysis_dir = args.output_root / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    population = pd.read_csv(run_dir / "population_summary.csv")
    endpoints = pd.read_csv(run_dir / "endpoint_counts.csv.gz")
    predicted = pd.read_csv(run_dir / "predicted_lox_frequency.csv")

    full = full_sample_table(population)
    rare = rarefaction(endpoints, population, args.rarefaction_reps, args.seed)
    comparison, rank_summary = diversity_comparison(full, rare)
    pareto = pareto_frontier(population, rare)
    backtest = backtesting(predicted)
    holdout_frame = holdout(predicted)
    null_summary, null_distribution = null_controls(predicted, backtest, args.null_reps, args.seed + 1)
    risk = gene_risk(run_dir, args.output_root)

    full.to_csv(analysis_dir / "full_sample_diversity.csv", index=False)
    rare.to_csv(analysis_dir / "rarefied_diversity.csv", index=False)
    comparison.to_csv(analysis_dir / "diversity_comparison.csv", index=False)
    rank_summary.to_csv(analysis_dir / "diversity_rank_consistency.csv", index=False)
    pareto.to_csv(analysis_dir / "survival_diversity_pareto_frontier.csv", index=False)
    backtest.to_csv(analysis_dir / "backtesting_report.csv", index=False)
    holdout_frame.to_csv(analysis_dir / "chromosome_wise_holdout.csv", index=False)
    null_summary.to_csv(analysis_dir / "null_control_summary.csv", index=False)
    null_distribution.to_csv(analysis_dir / "null_control_distributions.csv.gz", index=False, compression="gzip")
    risk.to_csv(analysis_dir / "exact_coordinate_gene_risk.csv", index=False)

    structural = rare[rare["endpoint_type"] == "structural"].sort_values(
        "effective_shannon_diversity_mean", ascending=False
    )
    best_structural = structural.iloc[0].to_dict()
    best_hic = backtest[backtest["model_key"] == "partial_hic_fallback"].sort_values(
        "pearson_r", ascending=False
    ).iloc[0].to_dict()
    summary = {
        "scenario_count": len(population),
        "common_rarefied_survivor_depth": int(rare["rarefied_survivor_depth"].min()),
        "rarefaction_replicates": args.rarefaction_reps,
        "best_structural_rarefied_diversity": best_structural,
        "best_partially_hic_informed_pearson": best_hic,
        "pareto_scenario_count": int(pareto["pareto_optimal"].sum()),
        "null_control_replicates": args.null_reps,
        "gene_risk_rows": len(risk),
    }
    (analysis_dir / "analysis_summary.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
