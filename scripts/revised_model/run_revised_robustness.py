from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run revised seed and essentiality-gate robustness scenarios.")
    parser.add_argument("--project-root", type=Path, default=REPOSITORY_ROOT)
    parser.add_argument("--output-root", type=Path, default=REPOSITORY_ROOT / "outputs" / "revised_model")
    parser.add_argument("--n-trajectories", type=int, default=10_000)
    parser.add_argument("--n-steps", type=int, default=20)
    parser.add_argument("--workers", type=int, default=max(1, min(8, (os.cpu_count() or 4) - 1)))
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def robustness_grid(output_root: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    seeds = [202606, 202607, 202608, 202609, 202610]
    for model_key in ["linear_distance", "partial_hic_fallback"]:
        for p_event in [0.09, 0.12]:
            for seed in seeds:
                rows.append(
                    {
                        "analysis_family": "five_seed_robustness",
                        "model_key": model_key,
                        "alpha": 2.0,
                        "p_event": p_event,
                        "gate_name": "strict",
                        "seed": seed,
                    }
                )
    gate_summary = pd.read_csv(output_root / "inputs" / "survival_gate_summary.csv")
    for gate_name in gate_summary["gate_name"].astype(str):
        for p_event in [0.09, 0.12]:
            rows.append(
                {
                    "analysis_family": "essentiality_gate_sensitivity",
                    "model_key": "linear_distance",
                    "alpha": 2.0,
                    "p_event": p_event,
                    "gate_name": gate_name,
                    "seed": 202606,
                }
            )
    for index, row in enumerate(rows, start=1):
        row["robustness_index"] = index
        row["run_subdir"] = f"robustness/scenarios/robustness_{index:03d}"
    return rows


def run_worker(args: argparse.Namespace, row: dict[str, object], log_dir: Path) -> dict[str, object]:
    script = Path(__file__).with_name("run_revised_simulation.py")
    index = int(row["robustness_index"])
    log_path = log_dir / f"robustness_{index:03d}.log"
    command = [
        sys.executable,
        str(script),
        "--project-root", str(args.project_root),
        "--output-root", str(args.output_root),
        "--n-trajectories", str(args.n_trajectories),
        "--n-steps", str(args.n_steps),
        "--seed", str(row["seed"]),
        "--model-key", str(row["model_key"]),
        "--alpha", str(row["alpha"]),
        "--p-event", str(row["p_event"]),
        "--gate-name", str(row["gate_name"]),
        "--run-subdir", str(row["run_subdir"]),
        "--reuse-gate-files",
    ]
    started = time.perf_counter()
    with log_path.open("w", encoding="utf-8") as log:
        completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True)
    return {
        **row,
        "return_code": completed.returncode,
        "elapsed_seconds": time.perf_counter() - started,
        "log_file": str(log_path),
        "command": json.dumps(command),
    }


def rarefy_robustness(output_root: Path, grid: list[dict[str, object]], depth: int, reps: int = 100) -> pd.DataFrame:
    rng = np.random.default_rng(20260718)
    rows: list[dict[str, object]] = []
    for row in grid:
        counts_path = output_root / str(row["run_subdir"]) / "endpoint_counts.csv.gz"
        counts_frame = pd.read_csv(counts_path)
        for endpoint_type, group in counts_frame.groupby("endpoint_type"):
            counts = group["count"].to_numpy(dtype=np.int64)
            values = []
            for _ in range(reps):
                sampled = rng.multivariate_hypergeometric(counts, depth)
                sampled = sampled[sampled > 0]
                probabilities = sampled / sampled.sum()
                shannon = float(-(probabilities * np.log(probabilities)).sum())
                values.append(
                    (
                        len(sampled),
                        math.exp(shannon),
                        float(1.0 / np.sum(probabilities**2)),
                    )
                )
            matrix = np.asarray(values, dtype=float)
            rows.append(
                {
                    "robustness_index": row["robustness_index"],
                    "analysis_family": row["analysis_family"],
                    "gate_name": row["gate_name"],
                    "model_key": row["model_key"],
                    "alpha": row["alpha"],
                    "p_event": row["p_event"],
                    "seed": row["seed"],
                    "endpoint_type": endpoint_type,
                    "rarefied_survivor_depth": depth,
                    "rarefaction_replicates": reps,
                    "unique_endpoints_mean": float(matrix[:, 0].mean()),
                    "effective_shannon_mean": float(matrix[:, 1].mean()),
                    "effective_shannon_ci_low": float(np.quantile(matrix[:, 1], 0.025)),
                    "effective_shannon_ci_high": float(np.quantile(matrix[:, 1], 0.975)),
                    "inverse_simpson_mean": float(matrix[:, 2].mean()),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    robustness_root = args.output_root / "robustness"
    log_dir = robustness_root / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    grid = robustness_grid(args.output_root)
    pd.DataFrame(grid).to_csv(robustness_root / "robustness_scenario_grid.csv", index=False)
    records: list[dict[str, object]] = []
    manifest_partial = robustness_root / "robustness_execution_manifest.partial.csv"
    successful: set[int] = set()
    if args.resume and manifest_partial.exists():
        previous = pd.read_csv(manifest_partial)
        previous = previous.sort_values("robustness_index").drop_duplicates("robustness_index", keep="last")
        previous = previous[previous["return_code"] == 0]
        records = previous.to_dict("records")
        successful = set(previous["robustness_index"].astype(int))
    pending = [row for row in grid if int(row["robustness_index"]) not in successful]
    print(f"successful_existing={len(successful)} pending={len(pending)} workers={args.workers}", flush=True)
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(run_worker, args, row, log_dir): row for row in pending}
        for future in as_completed(futures):
            record = future.result()
            records.append(record)
            print(
                f"robustness {record['robustness_index']}/{len(grid)} "
                f"family={record['analysis_family']} return={record['return_code']} "
                f"elapsed={record['elapsed_seconds']:.1f}s",
                flush=True,
            )
            pd.DataFrame(records).sort_values("robustness_index").to_csv(
                manifest_partial, index=False
            )
    manifest = pd.DataFrame(records).sort_values("robustness_index")
    manifest.to_csv(robustness_root / "robustness_execution_manifest.csv", index=False)
    failed = manifest[manifest["return_code"] != 0]
    if not failed.empty:
        raise RuntimeError(f"Robustness scenarios failed: {failed['robustness_index'].tolist()}")

    summaries: list[pd.DataFrame] = []
    for row in grid:
        path = args.output_root / str(row["run_subdir"]) / "population_summary.csv"
        frame = pd.read_csv(path)
        for key in ["analysis_family", "gate_name", "robustness_index"]:
            frame[key] = row[key]
        summaries.append(frame)
    combined = pd.concat(summaries, ignore_index=True)
    combined.to_csv(robustness_root / "robustness_population_summary.csv", index=False)

    common_depth = int(combined["gate_passing_trajectories"].min())
    rare = rarefy_robustness(args.output_root, grid, common_depth, reps=100)
    rare.to_csv(robustness_root / "robustness_rarefied_diversity.csv", index=False)
    structural_rare = rare[rare["endpoint_type"] == "structural"][
        ["robustness_index", "effective_shannon_mean", "effective_shannon_ci_low", "effective_shannon_ci_high"]
    ]
    combined = combined.merge(structural_rare, on="robustness_index", how="left")

    seed_summary = (
        combined[combined["analysis_family"] == "five_seed_robustness"]
        .groupby(["model_key", "alpha", "p_event"], as_index=False)
        .agg(
            replicate_count=("seed", "nunique"),
            gate_passing_fraction_mean=("essentiality_gate_passing_fraction", "mean"),
            gate_passing_fraction_sd=("essentiality_gate_passing_fraction", "std"),
            structural_effective_shannon_mean=("structural_effective_shannon_full_sample", "mean"),
            structural_effective_shannon_sd=("structural_effective_shannon_full_sample", "std"),
            rarefied_structural_effective_shannon_mean=("effective_shannon_mean", "mean"),
            rarefied_structural_effective_shannon_sd=("effective_shannon_mean", "std"),
        )
    )
    seed_summary.to_csv(robustness_root / "five_seed_robustness_summary.csv", index=False)
    gate_summary = combined[combined["analysis_family"] == "essentiality_gate_sensitivity"].copy()
    gate_summary.to_csv(robustness_root / "essentiality_gate_sensitivity_summary.csv", index=False)
    print(json.dumps({"scenario_count": len(grid), "all_return_codes_zero": True, "common_rarefied_depth": common_depth}, indent=2))


if __name__ == "__main__":
    main()
