from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

from run_revised_simulation import build_gate_sets, scenario_grid


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the 88-scenario revised model in isolated worker processes.")
    parser.add_argument("--project-root", type=Path, default=REPOSITORY_ROOT)
    parser.add_argument("--output-root", type=Path, default=REPOSITORY_ROOT / "outputs" / "revised_model")
    parser.add_argument("--n-trajectories", type=int, default=10_000)
    parser.add_argument("--n-steps", type=int, default=20)
    parser.add_argument("--seed", type=int, default=202606)
    parser.add_argument("--workers", type=int, default=max(1, min(8, (os.cpu_count() or 4) - 1)))
    parser.add_argument("--resume", action="store_true", help="Reuse successful entries and rerun only missing or failed scenarios.")
    return parser.parse_args()


def run_worker(args: argparse.Namespace, index: int, log_dir: Path) -> dict[str, object]:
    script = Path(__file__).with_name("run_revised_simulation.py")
    log_path = log_dir / f"scenario_{index:03d}.log"
    command = [
        sys.executable,
        str(script),
        "--project-root", str(args.project_root),
        "--output-root", str(args.output_root),
        "--n-trajectories", str(args.n_trajectories),
        "--n-steps", str(args.n_steps),
        "--seed", str(args.seed),
        "--scenario-index", str(index),
        "--reuse-gate-files",
    ]
    start = time.perf_counter()
    with log_path.open("w", encoding="utf-8") as log:
        completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, text=True)
    return {
        "scenario_index": index,
        "return_code": completed.returncode,
        "elapsed_seconds": time.perf_counter() - start,
        "log_file": str(log_path),
        "command": json.dumps(command),
    }


def merge_csv(paths: list[Path], destination: Path, compression: str | None = None) -> None:
    frames = [pd.read_csv(path, low_memory=False) for path in paths]
    pd.concat(frames, ignore_index=True).to_csv(destination, index=False, compression=compression)


def main() -> None:
    args = parse_args()
    inputs = args.output_root / "inputs"
    required = [
        inputs / "lox_sites_validated.csv",
        inputs / "reference_segments.csv",
        inputs / "hic_direct_pair_weights.csv",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Revised inputs must be built before the formal run: {missing}")
    build_gate_sets(args.output_root, args.project_root)
    run_root = args.output_root / "main_run"
    log_dir = run_root / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    scenarios = scenario_grid(False, False)
    records: list[dict[str, object]] = []
    manifest_partial = run_root / "scenario_execution_manifest.partial.csv"
    successful: set[int] = set()
    if args.resume and manifest_partial.exists():
        previous = pd.read_csv(manifest_partial)
        previous = previous.sort_values("scenario_index").drop_duplicates("scenario_index", keep="last")
        previous = previous[previous["return_code"] == 0]
        records = previous.to_dict("records")
        successful = set(previous["scenario_index"].astype(int))
    pending = [index for index in range(1, len(scenarios) + 1) if index not in successful]
    print(f"successful_existing={len(successful)} pending={len(pending)} workers={args.workers}", flush=True)
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(run_worker, args, index, log_dir): index
            for index in pending
        }
        for future in as_completed(futures):
            record = future.result()
            records.append(record)
            print(
                f"scenario {record['scenario_index']}/{len(scenarios)} return={record['return_code']} "
                f"elapsed={record['elapsed_seconds']:.1f}s",
                flush=True,
            )
            pd.DataFrame(records).sort_values("scenario_index").to_csv(
                manifest_partial, index=False
            )
    manifest = pd.DataFrame(records).sort_values("scenario_index")
    manifest.to_csv(run_root / "scenario_execution_manifest.csv", index=False)
    failed = manifest[manifest["return_code"] != 0]
    if not failed.empty:
        raise RuntimeError(f"Formal run failed for scenarios: {failed['scenario_index'].tolist()}")

    scenario_dirs = [run_root / "scenarios" / f"scenario_{index:03d}" for index in range(1, len(scenarios) + 1)]
    merge_csv([path / "population_summary.csv" for path in scenario_dirs], run_root / "population_summary.csv")
    merge_csv([path / "predicted_lox_frequency.csv" for path in scenario_dirs], run_root / "predicted_lox_frequency.csv")
    merge_csv(
        [path / "endpoint_counts.csv.gz" for path in scenario_dirs],
        run_root / "endpoint_counts.csv.gz",
        compression="gzip",
    )
    summary = {
        "scenario_count": len(scenarios),
        "workers": args.workers,
        "n_trajectories_per_scenario": args.n_trajectories,
        "all_return_codes_zero": True,
        "total_worker_seconds": float(manifest["elapsed_seconds"].sum()),
        "maximum_worker_seconds": float(manifest["elapsed_seconds"].max()),
    }
    (run_root / "parallel_run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
