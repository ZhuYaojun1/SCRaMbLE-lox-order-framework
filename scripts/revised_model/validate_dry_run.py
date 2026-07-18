from __future__ import annotations

import argparse
import gzip
import json
import sqlite3
import zlib
from pathlib import Path

import pandas as pd

from common import combination_count
from segment_engine import load_reference_genome


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate revised model inputs and dry-run outputs.")
    parser.add_argument("--output-root", type=Path, default=REPOSITORY_ROOT / "outputs" / "revised_model")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    inputs = args.output_root / "inputs"
    dry = args.output_root / "dry_run"
    checks: list[dict[str, object]] = []

    def check(name: str, passed: bool, value: object, expected: object, notes: str = "") -> None:
        checks.append(
            {
                "check": name,
                "status": "PASS" if passed else "FAIL",
                "value": value,
                "expected": expected,
                "notes": notes,
            }
        )

    lox = pd.read_csv(inputs / "lox_sites_validated.csv")
    counts = lox.groupby("chromosome").size().to_dict()
    check("validated_lox_count", len(lox) == 488, len(lox), 488)
    check(
        "chromosome_lox_counts",
        counts == {"SynII": 273, "SynIII": 100, "SynIXR": 44, "SynVI": 71},
        counts,
        {"SynII": 273, "SynIII": 100, "SynIXR": 44, "SynVI": 71},
    )
    pair_count = sum(combination_count(int(value)) for value in counts.values())
    check("candidate_pair_count", pair_count == 45_509, pair_count, 45_509)
    check("coordinate_resolved_lox", int(lox["position_bp"].notna().sum()) == 487, int(lox["position_bp"].notna().sum()), 487)

    audit = pd.read_csv(inputs / "lox_input_audit.csv")
    removed = int(audit["removed_blank_records"].sum())
    check("removed_padded_blank_records", removed == 604, removed, 604)
    predicted = pd.read_csv(dry / "predicted_lox_frequency.csv")
    rows_per_scenario = predicted.groupby(["model_key", "alpha", "p_event", "seed"]).size().unique().tolist()
    check("predicted_rows_per_scenario", rows_per_scenario == [488], rows_per_scenario, [488])

    hic = pd.read_csv(inputs / "hic_pair_coverage.csv")
    direct_non_synii = int(hic[(hic["chromosome"].isin(["SynIII", "SynVI", "SynIXR"]))]["direct_hic_pairs"].sum())
    check("non_synii_direct_hic_pairs", direct_non_synii == 0, direct_non_synii, 0)

    reference = load_reference_genome(inputs)
    reference.validate()
    segments = sum(len(topology.segments()) for topology in reference.chromosomes.values())
    active_lox = sum(len(topology.lox()) for topology in reference.chromosomes.values())
    check("reference_active_lox", active_lox == 488, active_lox, 488)
    check("reference_segment_count", segments == 492, segments, 492)

    population = pd.read_csv(dry / "population_summary.csv")
    trajectory = pd.read_csv(dry / "trajectory_summary.csv.gz")
    accepted = int(population["accepted_events"].sum())
    trajectory_events = int(trajectory["accepted_events"].sum())
    check("accepted_event_reconciliation", accepted == trajectory_events, trajectory_events, accepted)
    initialized = int(population["initialized_trajectories"].sum())
    check("trajectory_record_reconciliation", len(trajectory) == initialized, len(trajectory), initialized)
    terminal = int(population["gate_passing_trajectories"].sum() + population["gate_failing_trajectories"].sum())
    check("terminal_classification_reconciliation", terminal == initialized, terminal, initialized)

    endpoints = pd.read_csv(dry / "endpoint_counts.csv.gz")
    endpoint_totals = endpoints.groupby(["model_key", "alpha", "p_event", "seed", "endpoint_type"])["count"].sum()
    expected_totals = population.set_index(["model_key", "alpha", "p_event", "seed"])["gate_passing_trajectories"]
    endpoints_ok = all(
        int(value) == int(expected_totals.loc[key[:4]]) for key, value in endpoint_totals.items()
    )
    check("endpoint_count_reconciliation", endpoints_ok, "all endpoint types checked", "each equals gate-passing trajectories")

    connection = sqlite3.connect(dry / "endpoint_catalog.sqlite")
    sample = connection.execute("SELECT structural_hash, canonical_zlib FROM endpoint_structure LIMIT 20").fetchall()
    connection.close()
    reconstructable = all(
        zlib.decompress(blob).decode("utf-8").startswith("segment-topology-v1||") for _, blob in sample
    )
    check("endpoint_catalog_reconstructable", reconstructable, f"sampled={len(sample)}", "all canonical strings decompress")

    test_log = args.output_root / "unit_test_results.txt"
    check("event_unit_tests", test_log.exists() and "Ran 7 tests" in test_log.read_text(encoding="utf-8", errors="ignore") and "OK" in test_log.read_text(encoding="utf-8", errors="ignore"), str(test_log), "seven tests pass")

    failed = [row for row in checks if row["status"] == "FAIL"]
    table = pd.DataFrame(checks)
    table.to_csv(args.output_root / "dry_run_validation.csv", index=False)
    lines = [
        "# Dry-run validation",
        "",
        f"Overall status: **{'PASS' if not failed else 'FAIL'}**",
        "",
        "| Check | Status | Value | Expected | Notes |",
        "| --- | --- | --- | --- | --- |",
    ]
    for row in checks:
        lines.append(
            f"| {row['check']} | {row['status']} | {str(row['value']).replace('|', '/')} | "
            f"{str(row['expected']).replace('|', '/')} | {str(row['notes']).replace('|', '/')} |"
        )
    lines.extend(
        [
            "",
            "## Coordinate scope",
            "",
            "SynII, SynIII and SynVI have complete motif-supported lox coordinates for their Fig. 2C records. "
            "SynIXR has 44 non-empty Fig. 2C records but 43 verifiable local-FASTA loxPsym motifs. The unmatched "
            "record is retained only as an order-defined topology node; SynIXR is excluded from gene-level risk mapping.",
            "",
            "## Event semantics",
            "",
            "Deletion removes the current interval and one boundary copy; inversion reverses current segment order and "
            "orientation; duplication creates new segment and lox copy instances. Pair eligibility is rebuilt from the "
            "current topology at every accepted event.",
            "",
            "Formal execution is permitted only when every automated check above passes.",
        ]
    )
    (args.output_root / "dry_run_validation.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PASS" if not failed else "FAIL", "checks": len(checks), "failed": len(failed)}, indent=2))
    if failed:
        raise RuntimeError(f"Dry-run validation failed: {[row['check'] for row in failed]}")


if __name__ == "__main__":
    main()
