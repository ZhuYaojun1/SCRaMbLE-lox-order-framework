from __future__ import annotations

import argparse
import json
import re
import sqlite3
import zlib
from pathlib import Path

import pandas as pd


EVENT_COPY = re.compile(r"@e\d+")
REFERENCE_CHROMOSOME = re.compile(r"(?:^|\|\|)([^|]+)\|REFERENCE:")
REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Classify provenance hashes merged by unlabeled topology.")
    parser.add_argument("--output-root", type=Path, default=REPOSITORY_ROOT / "outputs" / "revised_model")
    return parser.parse_args()


def load_payloads(database: Path, hashes: list[str]) -> dict[str, str]:
    connection = sqlite3.connect(database)
    try:
        placeholders = ",".join("?" for _ in hashes)
        rows = connection.execute(
            f"SELECT structural_hash, canonical_zlib FROM endpoint_structure WHERE structural_hash IN ({placeholders})",
            hashes,
        ).fetchall()
        return {endpoint_hash: zlib.decompress(blob).decode("utf-8") for endpoint_hash, blob in rows}
    finally:
        connection.close()


def main() -> None:
    args = parse_args()
    analysis = args.output_root / "analysis"
    mapping = pd.read_csv(analysis / "endpoint_provenance_to_unlabeled_topology.csv.gz")
    rows = []
    for (scenario_index, topology_hash), group in mapping.groupby(
        ["scenario_index", "unlabeled_topology_hash"], sort=True
    ):
        if len(group) <= 1:
            continue
        hashes = group["provenance_endpoint_hash"].astype(str).tolist()
        scenario_dir = args.output_root / "main_run" / "scenarios" / f"scenario_{int(scenario_index):03d}"
        payloads = load_payloads(scenario_dir / "endpoint_catalog.sqlite", hashes)
        has_event_copy = any(EVENT_COPY.search(payload) for payload in payloads.values())
        reference_sets = [frozenset(REFERENCE_CHROMOSOME.findall(payload)) for payload in payloads.values()]
        if has_event_copy:
            classification = "event_copy_identity_alias"
        elif len(set(reference_sets)) > 1:
            classification = "reference_token_vs_expanded_payload_alias"
        else:
            classification = "other_canonical_representation_alias"
        first = group.iloc[0]
        rows.append(
            {
                "scenario_index": int(scenario_index),
                "model_key": first["model_key"],
                "alpha": first["alpha"],
                "p_event": first["p_event"],
                "unlabeled_topology_hash": topology_hash,
                "classification": classification,
                "provenance_hash_count": len(hashes),
                "excess_provenance_hashes": len(hashes) - 1,
                "trajectory_count": int(group["trajectory_count"].sum()),
                "contains_event_copy_id": has_event_copy,
                "reference_token_chromosome_sets": ";".join(",".join(sorted(value)) for value in reference_sets),
                "provenance_hashes": ";".join(hashes),
            }
        )
    aliases = pd.DataFrame(rows)
    aliases.to_csv(analysis / "endpoint_alias_classification.csv", index=False)
    summary = (
        aliases.groupby("classification", as_index=False)
        .agg(
            topology_groups=("unlabeled_topology_hash", "count"),
            excess_provenance_hashes=("excess_provenance_hashes", "sum"),
            trajectories_in_groups=("trajectory_count", "sum"),
        )
        .sort_values("topology_groups", ascending=False)
    )
    summary.to_csv(analysis / "endpoint_alias_classification_summary.csv", index=False)
    print(json.dumps({"merged_topology_groups": len(aliases), "classification": summary.to_dict("records")}, indent=2))


if __name__ == "__main__":
    main()
