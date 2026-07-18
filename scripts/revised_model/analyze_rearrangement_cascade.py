from __future__ import annotations

import argparse
import csv
import heapq
import json
import lzma
import math
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from run_revised_simulation import event_probabilities, load_gate_sets
from segment_engine import GenomeState, PairSampler, clone_reference, load_reference_genome


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
PRECURSOR_TYPES = {"inversion", "duplication"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Replay Gate-failing trajectories and quantify rearrangement cascade effects."
    )
    parser.add_argument(
        "--event-root",
        type=Path,
        required=True,
        help="Directory containing scenario_*/event_log.csv.gz from the formal run.",
    )
    parser.add_argument(
        "--inputs-dir",
        type=Path,
        default=REPOSITORY_ROOT / "outputs" / "revised_model" / "inputs",
    )
    parser.add_argument("--project-root", type=Path, default=REPOSITORY_ROOT)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPOSITORY_ROOT / "outputs" / "revised_model" / "cascade_analysis",
    )
    parser.add_argument("--gate-name", default="strict")
    parser.add_argument(
        "--total-hazard-sample-per-type",
        type=int,
        default=50,
        help="Per-scenario sample of each precursor type for exact all-pair lethal-hazard calculation.",
    )
    parser.add_argument("--seed", type=int, default=20260718)
    parser.add_argument("--scenario-limit", type=int, default=None)
    return parser.parse_args()


def split_ids(value: Any) -> list[str]:
    if value is None:
        return []
    text = str(value)
    if not text or text == "nan":
        return []
    return [token for token in text.split(";") if token]


def classify_history(events: list[str]) -> str:
    if not events:
        return "deletion_only"
    kinds = set(events)
    has_inversion = "inversion" in kinds
    has_duplication = "duplication" in kinds
    has_deletion = "deletion" in kinds
    if has_deletion and not (has_inversion or has_duplication):
        return "earlier_deletion_only"
    if has_deletion and (has_inversion or has_duplication):
        return "earlier_deletion_with_nondeletion"
    if has_inversion and has_duplication:
        return "mixed_inversion_duplication"
    if has_inversion:
        return "inversion_only_history"
    if has_duplication:
        return "duplication_only_history"
    return "other"


def wilson_interval(successes: int, total: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if total <= 0:
        return float("nan"), float("nan")
    proportion = successes / total
    denominator = 1.0 + z * z / total
    centre = (proportion + z * z / (2.0 * total)) / denominator
    half = z * math.sqrt(proportion * (1.0 - proportion) / total + z * z / (4.0 * total * total)) / denominator
    return centre - half, centre + half


class CascadeProbability:
    def __init__(self, direct: pd.DataFrame) -> None:
        synii_ids = sorted(
            set(direct["lox_id_1"].astype(str)) | set(direct["lox_id_2"].astype(str)),
            key=lambda value: int(value.rsplit("_", 1)[-1]),
        )
        self.reference_ids = [f"{value}@0" for value in synii_ids]
        self.reference_index = {value: index for index, value in enumerate(self.reference_ids)}
        size = len(self.reference_ids)
        self.direct_matrix = np.zeros((size, size), dtype=np.float64)
        edge_u: list[int] = []
        edge_v: list[int] = []
        edge_w: list[float] = []
        for row in direct.to_dict("records"):
            first = self.reference_index[f"{row['lox_id_1']}@0"]
            second = self.reference_index[f"{row['lox_id_2']}@0"]
            weight = float(row["direct_contact_weight"])
            self.direct_matrix[first, second] = weight
            self.direct_matrix[second, first] = weight
            edge_u.append(first)
            edge_v.append(second)
            edge_w.append(weight)
        self.edge_u = np.asarray(edge_u, dtype=np.int16)
        self.edge_v = np.asarray(edge_v, dtype=np.int16)
        self.edge_w = np.asarray(edge_w, dtype=np.float64)
        self.direct_total_cache: dict[tuple[int, ...], float] = {}

    def active_reference_positions(self, state: GenomeState) -> tuple[np.ndarray, tuple[int, ...]]:
        positions = np.full(len(self.reference_ids), -1, dtype=np.int16)
        active: list[int] = []
        for position, lox in enumerate(state.chromosomes["SynII"].lox()):
            reference_index = self.reference_index.get(lox.copy_id)
            if reference_index is not None:
                positions[reference_index] = position
                active.append(reference_index)
        return positions, tuple(sorted(active))

    def eligible_direct_total(self, state: GenomeState) -> float:
        _, active = self.active_reference_positions(state)
        cached = self.direct_total_cache.get(active)
        if cached is not None:
            return cached
        indices = np.asarray(active, dtype=int)
        if indices.size < 2:
            total = 0.0
        else:
            total = float(self.direct_matrix[np.ix_(indices, indices)].sum() / 2.0)
        self.direct_total_cache[active] = total
        return total

    def pair_metrics(
        self,
        state: GenomeState,
        sampler: PairSampler,
        gate_counts: Counter[str],
        gate_genes: set[str],
        chromosome: str,
        left_copy_id: str,
        right_copy_id: str,
        p_event: float,
        deletion_probability: float,
    ) -> dict[str, Any]:
        topology = state.chromosomes[chromosome]
        indices = topology.lox_index_by_copy_id()
        if left_copy_id not in indices or right_copy_id not in indices:
            return {
                "exists": False,
                "distance": np.nan,
                "proposal_probability": 0.0,
                "deletion_is_gate_failing": False,
                "lethal_deletion_hazard": 0.0,
                "spanned_gate_genes": "",
            }
        left, right = sorted((indices[left_copy_id], indices[right_copy_id]))
        if left == right:
            raise RuntimeError("A lox pair collapsed to one active copy.")
        distance = right - left
        _, _, fallback_total = sampler._fallback_totals(state)
        if sampler.model_name == "uniform_random":
            fallback_weight = 1.0
            denominator = fallback_total
            direct_weight = 0.0
        else:
            fallback_weight = 1.0 / ((1.0 + distance) ** sampler.alpha)
            if sampler.model_name == "partial_hic_fallback":
                key = (chromosome, *sorted((left_copy_id, right_copy_id)))
                direct_weight = float(sampler.direct_lookup.get(key, 0.0))
                denominator = fallback_total + self.eligible_direct_total(state)
            else:
                direct_weight = 0.0
                denominator = fallback_total
        probability = (fallback_weight + direct_weight) / denominator if denominator > 0 else 0.0
        removed_counts: Counter[str] = Counter()
        for segment in topology.segments()[left + 1 : right + 1]:
            removed_counts.update(gene for gene in segment.orf_ids if gene in gate_genes)
        spanned = sorted(
            gene
            for gene, removed in removed_counts.items()
            if gate_counts.get(gene, 0) > 0 and removed >= gate_counts[gene]
        )
        lethal = bool(spanned)
        return {
            "exists": True,
            "distance": distance,
            "proposal_probability": probability,
            "deletion_is_gate_failing": lethal,
            "lethal_deletion_hazard": p_event * deletion_probability * probability if lethal else 0.0,
            "spanned_gate_genes": ";".join(spanned),
        }

    def total_lethal_hazard(
        self,
        state: GenomeState,
        sampler: PairSampler,
        gate_genes: set[str],
        p_event: float,
        deletion_probability: float,
    ) -> dict[str, float]:
        fallback_lethal_mass = 0.0
        direct_lethal_mass = 0.0
        synii_thresholds: np.ndarray | None = None
        for chromosome, topology in state.chromosomes.items():
            segments = topology.segments()
            n_lox = len(topology.lox())
            if n_lox < 2:
                continue
            extents: dict[str, list[int]] = {}
            for segment_index, segment in enumerate(segments):
                for gene in segment.orf_ids:
                    if gene not in gate_genes:
                        continue
                    if gene not in extents:
                        extents[gene] = [segment_index, segment_index]
                    else:
                        extents[gene][0] = min(extents[gene][0], segment_index)
                        extents[gene][1] = max(extents[gene][1], segment_index)
            best_at_start = np.full(n_lox + 1, n_lox + 1, dtype=np.int32)
            for minimum, maximum in extents.values():
                if minimum > 0:
                    best_at_start[minimum] = min(best_at_start[minimum], maximum)
            suffix = np.minimum.accumulate(best_at_start[::-1])[::-1]
            thresholds = np.full(n_lox, n_lox + 1, dtype=np.int32)
            for left in range(n_lox - 1):
                threshold = int(suffix[left + 1])
                thresholds[left] = threshold
                first_distance = max(1, threshold - left)
                last_distance = n_lox - 1 - left
                if first_distance > last_distance:
                    continue
                if sampler.model_name == "uniform_random":
                    fallback_lethal_mass += last_distance - first_distance + 1
                else:
                    distances = np.arange(first_distance, last_distance + 1, dtype=float)
                    fallback_lethal_mass += float(np.sum((1.0 + distances) ** (-sampler.alpha)))
            if chromosome == "SynII":
                synii_thresholds = thresholds

        _, _, fallback_total = sampler._fallback_totals(state)
        direct_eligible_total = 0.0
        if sampler.model_name == "partial_hic_fallback":
            positions, active = self.active_reference_positions(state)
            active_u = positions[self.edge_u]
            active_v = positions[self.edge_v]
            eligible = (active_u >= 0) & (active_v >= 0)
            direct_eligible_total = float(self.edge_w[eligible].sum())
            if synii_thresholds is not None and np.any(eligible):
                left = np.minimum(active_u[eligible], active_v[eligible]).astype(int)
                right = np.maximum(active_u[eligible], active_v[eligible]).astype(int)
                lethal = right >= synii_thresholds[left]
                direct_lethal_mass = float(self.edge_w[eligible][lethal].sum())
        denominator = fallback_total + direct_eligible_total
        lethal_probability = (
            (fallback_lethal_mass + direct_lethal_mass) / denominator if denominator > 0 else 0.0
        )
        return {
            "total_lethal_pair_probability": lethal_probability,
            "total_lethal_deletion_hazard": p_event * deletion_probability * lethal_probability,
            "fallback_lethal_mass": fallback_lethal_mass,
            "direct_lethal_mass": direct_lethal_mass,
            "eligible_pair_mass": denominator,
        }


def interval_gate_counts(
    state: GenomeState,
    chromosome: str,
    left: int,
    right: int,
    gate_genes: set[str],
) -> Counter[str]:
    counts: Counter[str] = Counter()
    for segment in state.chromosomes[chromosome].segments()[left + 1 : right + 1]:
        counts.update(gene for gene in segment.orf_ids if gene in gate_genes)
    return counts


def finite_log2_ratio(after: float, before: float) -> float:
    if after > 0 and before > 0:
        return math.log2(after / before)
    return float("nan")


def mechanism(before: dict[str, Any], after: dict[str, Any]) -> str:
    if not before["exists"] and after["exists"]:
        return "final_pair_activation"
    if before["deletion_is_gate_failing"] and not after["deletion_is_gate_failing"]:
        return "essential_copy_buffering"
    if not before["deletion_is_gate_failing"] and after["deletion_is_gate_failing"]:
        return "lethal_span_activation"
    if before["exists"] and after["exists"]:
        distance_change = after["distance"] - before["distance"]
        probability_change = after["proposal_probability"] - before["proposal_probability"]
        if distance_change < 0 and probability_change > 0:
            return "distance_shortening"
        if distance_change > 0 and probability_change < 0:
            return "distance_lengthening"
        if distance_change == 0 and not math.isclose(probability_change, 0.0, abs_tol=1e-18):
            return "proposal_space_renormalization"
    return "no_direct_effect_on_final_pair"


def choose_total_hazard_sample(
    failed_events: pd.DataFrame,
    terminal_ids: pd.Series,
    sample_per_type: int,
    seed: int,
) -> set[tuple[int, int]]:
    candidates = failed_events.copy()
    candidates["terminal_event_id"] = candidates["trajectory_id"].map(terminal_ids)
    candidates = candidates[
        candidates["event_type"].isin(PRECURSOR_TYPES)
        & (candidates["event_id"] < candidates["terminal_event_id"])
    ]
    selected: set[tuple[int, int]] = set()
    rng = np.random.default_rng(seed)
    for event_type, group in candidates.groupby("event_type", sort=True):
        size = min(sample_per_type, len(group))
        if size <= 0:
            continue
        indices = rng.choice(len(group), size=size, replace=False)
        for row in group.iloc[indices].itertuples(index=False):
            selected.add((int(row.trajectory_id), int(row.event_id)))
    return selected


def summarize_precursors(frame: pd.DataFrame) -> pd.DataFrame:
    group_columns = ["model_key", "alpha", "p_event", "precursor_event_type"]
    rows: list[dict[str, Any]] = []
    grouped: list[tuple[Any, pd.DataFrame]] = list(frame.groupby(group_columns, sort=True))
    grouped.extend(
        [("ALL", group) for _, group in frame.groupby("precursor_event_type", sort=True)]
    )
    for key, group in grouped:
        if key == "ALL":
            continue
    for key, group in frame.groupby(group_columns, sort=True):
        rows.append(_precursor_summary_row(key, group, "scenario"))
    for event_type, group in frame.groupby("precursor_event_type", sort=True):
        rows.append(_precursor_summary_row(("ALL", np.nan, np.nan, event_type), group, "overall"))
    return pd.DataFrame(rows)


def _precursor_summary_row(key: tuple[Any, ...], group: pd.DataFrame, scope: str) -> dict[str, Any]:
    model_key, alpha, p_event, event_type = key
    n = len(group)
    probability_increased = int(group["proposal_probability_increased"].sum())
    hazard_increased = int(group["final_pair_lethal_hazard_increased"].sum())
    pair_activated = int(group["final_pair_activated"].sum())
    low, high = wilson_interval(hazard_increased, n)
    finite_ratio = pd.to_numeric(group["proposal_log2_ratio_nonzero"], errors="coerce").dropna()
    distances = pd.to_numeric(group["distance_change"], errors="coerce").dropna()
    sampled = group[group["total_hazard_sampled"]]
    total_ratio = pd.to_numeric(sampled["total_hazard_log2_ratio_nonzero"], errors="coerce").dropna()
    return {
        "scope": scope,
        "model_key": model_key,
        "alpha": alpha,
        "p_event": p_event,
        "precursor_event_type": event_type,
        "n_precursor_events": n,
        "n_trajectories": group["trajectory_key"].nunique(),
        "same_chromosome_fraction": float(group["same_chromosome_as_terminal"].mean()),
        "final_pair_activation_fraction": pair_activated / n,
        "proposal_probability_increase_fraction": probability_increased / n,
        "final_pair_lethal_hazard_increase_fraction": hazard_increased / n,
        "hazard_increase_ci_low_event_level": low,
        "hazard_increase_ci_high_event_level": high,
        "median_proposal_log2_ratio_when_nonzero": float(finite_ratio.median()) if len(finite_ratio) else np.nan,
        "median_distance_change_when_observable": float(distances.median()) if len(distances) else np.nan,
        "n_total_hazard_sampled": len(sampled),
        "total_lethal_hazard_increase_fraction_sample": float(sampled["total_lethal_hazard_increased"].mean()) if len(sampled) else np.nan,
        "median_total_hazard_log2_ratio_nonzero_sample": float(total_ratio.median()) if len(total_ratio) else np.nan,
    }


def summarize_trajectory_effects(frame: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for event_type in ["inversion", "duplication"]:
        count_column = f"{event_type}_precursor_count"
        selected = frame[frame[count_column] > 0]
        for model_key, group in [("ALL", selected), *list(selected.groupby("model_key", sort=True))]:
            n = len(group)
            successes = int(group[f"any_{event_type}_final_pair_hazard_increase"].sum())
            low, high = wilson_interval(successes, n)
            ratios = pd.to_numeric(group[f"max_{event_type}_proposal_log2_ratio_nonzero"], errors="coerce").dropna()
            rows.append(
                {
                    "model_key": model_key,
                    "precursor_event_type": event_type,
                    "n_trajectories": n,
                    "trajectory_fraction_with_any_final_pair_hazard_increase": successes / n if n else np.nan,
                    "trajectory_fraction_ci_low": low,
                    "trajectory_fraction_ci_high": high,
                    "trajectory_fraction_with_pair_activation": float(group[f"any_{event_type}_pair_activation"].mean()) if n else np.nan,
                    "median_of_trajectory_max_log2_probability_ratio": float(ratios.median()) if len(ratios) else np.nan,
                }
            )
    return pd.DataFrame(rows)


def markdown_table(frame: pd.DataFrame, columns: list[str]) -> str:
    selected = frame[columns].copy()
    lines = ["| " + " | ".join(columns) + " |", "|" + "|".join(["---"] * len(columns)) + "|"]
    for row in selected.to_dict("records"):
        values = []
        for column in columns:
            value = row[column]
            if isinstance(value, float):
                values.append("NA" if not math.isfinite(value) else f"{value:.4f}")
            else:
                values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    event_paths = sorted(args.event_root.glob("scenario_*/event_log.csv.gz"))
    if args.scenario_limit is not None:
        event_paths = event_paths[: args.scenario_limit]
    if not event_paths:
        raise FileNotFoundError(f"No scenario event logs found under {args.event_root}")

    reference = load_reference_genome(args.inputs_dir)
    gate_sets = load_gate_sets(args.inputs_dir.parent)
    gate_genes = gate_sets[args.gate_name]
    initial_counts = reference.orf_copy_numbers()
    if any(initial_counts.get(gene, 0) <= 0 for gene in gate_genes):
        raise RuntimeError("The selected gate contains genes absent from the reference topology.")
    event_probs = event_probabilities(args.project_root)
    deletion_probability = float(event_probs["deletion"])
    direct = pd.read_csv(args.inputs_dir / "hic_direct_pair_weights.csv")
    probability = CascadeProbability(direct)

    precursor_fields = [
        "scenario", "model_key", "alpha", "p_event", "seed", "trajectory_id", "trajectory_key",
        "precursor_event_id", "precursor_step", "precursor_event_type", "precursor_chromosome",
        "precursor_left_lox_copy_id", "precursor_right_lox_copy_id", "terminal_event_id",
        "terminal_step", "terminal_chromosome", "terminal_left_lox_copy_id",
        "terminal_right_lox_copy_id", "terminal_missing_gate_genes", "events_until_failure",
        "same_chromosome_as_terminal", "immediate_predecessor", "final_pair_exists_before",
        "final_pair_exists_after", "final_pair_activated", "final_pair_distance_before",
        "final_pair_distance_after", "distance_change", "final_pair_lethal_before",
        "final_pair_lethal_after", "final_pair_lethality_activated", "final_pair_buffered_after",
        "proposal_probability_before", "proposal_probability_after", "proposal_probability_change",
        "proposal_probability_increased", "proposal_log2_ratio_nonzero",
        "final_pair_lethal_hazard_before", "final_pair_lethal_hazard_after",
        "final_pair_lethal_hazard_increased", "cascade_mechanism", "total_hazard_sampled",
        "total_lethal_pair_probability_before", "total_lethal_pair_probability_after",
        "total_lethal_hazard_before", "total_lethal_hazard_after",
        "total_lethal_hazard_increased", "total_hazard_log2_ratio_nonzero",
    ]
    trajectory_fields = [
        "scenario", "model_key", "alpha", "p_event", "seed", "trajectory_id", "trajectory_key",
        "history_class", "events_before_terminal_deletion", "prior_deletion_count",
        "inversion_precursor_count", "duplication_precursor_count", "terminal_chromosome",
        "terminal_pair_weight_source", "terminal_pair_distance", "terminal_pair_proposal_probability",
        "terminal_missing_gate_genes", "any_inversion_pair_activation", "any_duplication_pair_activation",
        "any_inversion_final_pair_hazard_increase", "any_duplication_final_pair_hazard_increase",
        "max_inversion_proposal_log2_ratio_nonzero", "max_duplication_proposal_log2_ratio_nonzero",
    ]
    precursor_path = args.output_dir / "cascade_precursor_events.csv.xz"
    trajectory_path = args.output_dir / "cascade_trajectory_effects.csv.xz"
    precursor_handle = lzma.open(precursor_path, "wt", encoding="utf-8", newline="")
    trajectory_handle = lzma.open(trajectory_path, "wt", encoding="utf-8", newline="")
    precursor_writer = csv.DictWriter(precursor_handle, fieldnames=precursor_fields)
    trajectory_writer = csv.DictWriter(trajectory_handle, fieldnames=trajectory_fields)
    precursor_writer.writeheader()
    trajectory_writer.writeheader()

    sequence_counts: Counter[tuple[str, float, float, str]] = Counter()
    sequence_overall: Counter[str] = Counter()
    replay_audit: list[dict[str, Any]] = []
    top_events: list[tuple[float, int, dict[str, Any]]] = []
    top_counter = 0
    total_events = 0
    total_failures = 0
    total_precursors = 0
    total_replayed = 0
    start_all = time.perf_counter()

    try:
        for scenario_number, event_path in enumerate(event_paths, start=1):
            scenario_start = time.perf_counter()
            scenario = event_path.parent.name
            frame = pd.read_csv(event_path, keep_default_na=False)
            total_events += len(frame)
            terminal = frame[frame["gate_status_after_event"].eq("fail")].copy()
            if terminal["trajectory_id"].duplicated().any():
                raise RuntimeError(f"{scenario}: more than one terminal Gate-failing row per trajectory.")
            if not terminal["event_type"].eq("deletion").all():
                raise RuntimeError(f"{scenario}: a Gate-failing trajectory did not terminate with deletion.")
            terminal_ids = terminal.set_index("trajectory_id")["event_id"]
            failed_ids = set(terminal_ids.index.astype(int))
            failed = frame[frame["trajectory_id"].isin(failed_ids)].copy()
            failed.sort_values(["trajectory_id", "event_id"], inplace=True)
            total_failures += len(failed_ids)
            sample_keys = choose_total_hazard_sample(
                failed,
                terminal_ids,
                args.total_hazard_sample_per_type,
                args.seed + scenario_number,
            )
            model_key = str(frame["model_key"].iloc[0])
            alpha = float(frame["alpha"].iloc[0])
            p_event = float(frame["p_event"].iloc[0])
            seed = int(frame["seed"].iloc[0])
            sampler = PairSampler(
                model_key,
                alpha,
                direct if model_key == "partial_hic_fallback" else None,
            )
            scenario_precursors = 0
            scenario_replayed = 0
            scenario_mismatches = 0

            for trajectory_id, group in failed.groupby("trajectory_id", sort=False):
                rows = list(group.itertuples(index=False))
                terminal_row = rows[-1]
                prior_types = [str(row.event_type) for row in rows[:-1]]
                history = classify_history(prior_types)
                sequence_counts[(model_key, alpha, p_event, history)] += 1
                sequence_overall[history] += 1
                precursor_rows = [row for row in rows[:-1] if str(row.event_type) in PRECURSOR_TYPES]
                if not precursor_rows:
                    continue
                scenario_replayed += 1
                total_replayed += 1
                state = clone_reference(reference)
                gate_counts: Counter[str] = Counter(
                    {gene: int(initial_counts.get(gene, 0)) for gene in gate_genes}
                )
                event_effects: dict[str, list[dict[str, Any]]] = defaultdict(list)
                terminal_metrics: dict[str, Any] | None = None
                for row in rows:
                    chromosome = str(row.chromosome)
                    topology = state.chromosomes[chromosome]
                    index_lookup = topology.lox_index_by_copy_id()
                    left_id = str(row.left_lox_copy_id)
                    right_id = str(row.right_lox_copy_id)
                    if left_id not in index_lookup or right_id not in index_lookup:
                        raise RuntimeError(
                            f"{scenario} trajectory {trajectory_id}: replay could not locate event pair {left_id}, {right_id}."
                        )
                    left, right = sorted((index_lookup[left_id], index_lookup[right_id]))
                    affected_counts = interval_gate_counts(state, chromosome, left, right, gate_genes)
                    is_precursor = int(row.event_id) < int(terminal_row.event_id) and str(row.event_type) in PRECURSOR_TYPES
                    sampled = (int(trajectory_id), int(row.event_id)) in sample_keys
                    if is_precursor:
                        before = probability.pair_metrics(
                            state,
                            sampler,
                            gate_counts,
                            gate_genes,
                            str(terminal_row.chromosome),
                            str(terminal_row.left_lox_copy_id),
                            str(terminal_row.right_lox_copy_id),
                            p_event,
                            deletion_probability,
                        )
                        total_before = probability.total_lethal_hazard(
                            state, sampler, gate_genes, p_event, deletion_probability
                        ) if sampled else None
                    if int(row.event_id) == int(terminal_row.event_id):
                        terminal_metrics = probability.pair_metrics(
                            state,
                            sampler,
                            gate_counts,
                            gate_genes,
                            chromosome,
                            left_id,
                            right_id,
                            p_event,
                            deletion_probability,
                        )
                    outcome = state.apply_event(str(row.event_type), chromosome, left, right)
                    if outcome.left_lox_copy_id != left_id or outcome.right_lox_copy_id != right_id:
                        scenario_mismatches += 1
                    if str(row.event_type) == "deletion":
                        gate_counts.subtract(affected_counts)
                    elif str(row.event_type) == "duplication":
                        gate_counts.update(affected_counts)
                    if is_precursor:
                        after = probability.pair_metrics(
                            state,
                            sampler,
                            gate_counts,
                            gate_genes,
                            str(terminal_row.chromosome),
                            str(terminal_row.left_lox_copy_id),
                            str(terminal_row.right_lox_copy_id),
                            p_event,
                            deletion_probability,
                        )
                        total_after = probability.total_lethal_hazard(
                            state, sampler, gate_genes, p_event, deletion_probability
                        ) if sampled else None
                        probability_change = after["proposal_probability"] - before["proposal_probability"]
                        hazard_change = after["lethal_deletion_hazard"] - before["lethal_deletion_hazard"]
                        distance_change = (
                            after["distance"] - before["distance"]
                            if before["exists"] and after["exists"]
                            else np.nan
                        )
                        record = {
                            "scenario": scenario,
                            "model_key": model_key,
                            "alpha": alpha,
                            "p_event": p_event,
                            "seed": seed,
                            "trajectory_id": int(trajectory_id),
                            "trajectory_key": f"{scenario}:{int(trajectory_id)}",
                            "precursor_event_id": int(row.event_id),
                            "precursor_step": int(row.step),
                            "precursor_event_type": str(row.event_type),
                            "precursor_chromosome": chromosome,
                            "precursor_left_lox_copy_id": left_id,
                            "precursor_right_lox_copy_id": right_id,
                            "terminal_event_id": int(terminal_row.event_id),
                            "terminal_step": int(terminal_row.step),
                            "terminal_chromosome": str(terminal_row.chromosome),
                            "terminal_left_lox_copy_id": str(terminal_row.left_lox_copy_id),
                            "terminal_right_lox_copy_id": str(terminal_row.right_lox_copy_id),
                            "terminal_missing_gate_genes": str(terminal_row.missing_gate_genes),
                            "events_until_failure": int(terminal_row.event_id) - int(row.event_id),
                            "same_chromosome_as_terminal": chromosome == str(terminal_row.chromosome),
                            "immediate_predecessor": int(row.event_id) + 1 == int(terminal_row.event_id),
                            "final_pair_exists_before": before["exists"],
                            "final_pair_exists_after": after["exists"],
                            "final_pair_activated": (not before["exists"]) and after["exists"],
                            "final_pair_distance_before": before["distance"],
                            "final_pair_distance_after": after["distance"],
                            "distance_change": distance_change,
                            "final_pair_lethal_before": before["deletion_is_gate_failing"],
                            "final_pair_lethal_after": after["deletion_is_gate_failing"],
                            "final_pair_lethality_activated": (not before["deletion_is_gate_failing"]) and after["deletion_is_gate_failing"],
                            "final_pair_buffered_after": before["deletion_is_gate_failing"] and (not after["deletion_is_gate_failing"]),
                            "proposal_probability_before": before["proposal_probability"],
                            "proposal_probability_after": after["proposal_probability"],
                            "proposal_probability_change": probability_change,
                            "proposal_probability_increased": probability_change > 1e-18,
                            "proposal_log2_ratio_nonzero": finite_log2_ratio(after["proposal_probability"], before["proposal_probability"]),
                            "final_pair_lethal_hazard_before": before["lethal_deletion_hazard"],
                            "final_pair_lethal_hazard_after": after["lethal_deletion_hazard"],
                            "final_pair_lethal_hazard_increased": hazard_change > 1e-18,
                            "cascade_mechanism": mechanism(before, after),
                            "total_hazard_sampled": sampled,
                            "total_lethal_pair_probability_before": total_before["total_lethal_pair_probability"] if total_before else np.nan,
                            "total_lethal_pair_probability_after": total_after["total_lethal_pair_probability"] if total_after else np.nan,
                            "total_lethal_hazard_before": total_before["total_lethal_deletion_hazard"] if total_before else np.nan,
                            "total_lethal_hazard_after": total_after["total_lethal_deletion_hazard"] if total_after else np.nan,
                            "total_lethal_hazard_increased": (total_after["total_lethal_deletion_hazard"] - total_before["total_lethal_deletion_hazard"] > 1e-18) if total_before else False,
                            "total_hazard_log2_ratio_nonzero": finite_log2_ratio(total_after["total_lethal_deletion_hazard"], total_before["total_lethal_deletion_hazard"]) if total_before else np.nan,
                        }
                        precursor_writer.writerow(record)
                        event_effects[str(row.event_type)].append(record)
                        scenario_precursors += 1
                        total_precursors += 1
                        score = float(record["proposal_log2_ratio_nonzero"])
                        if record["final_pair_activated"]:
                            score = 1_000.0 + math.log10(max(record["proposal_probability_after"], 1e-300))
                        if math.isfinite(score):
                            top_counter += 1
                            item = (score, top_counter, record)
                            if len(top_events) < 200:
                                heapq.heappush(top_events, item)
                            elif score > top_events[0][0]:
                                heapq.heapreplace(top_events, item)

                missing_calculated = sorted(gene for gene in gate_genes if gate_counts.get(gene, 0) <= 0)
                missing_logged = sorted(split_ids(terminal_row.missing_gate_genes))
                if missing_calculated != missing_logged:
                    scenario_mismatches += 1
                if terminal_metrics is None or not terminal_metrics["deletion_is_gate_failing"]:
                    scenario_mismatches += 1
                trajectory_record: dict[str, Any] = {
                    "scenario": scenario,
                    "model_key": model_key,
                    "alpha": alpha,
                    "p_event": p_event,
                    "seed": seed,
                    "trajectory_id": int(trajectory_id),
                    "trajectory_key": f"{scenario}:{int(trajectory_id)}",
                    "history_class": history,
                    "events_before_terminal_deletion": len(prior_types),
                    "prior_deletion_count": prior_types.count("deletion"),
                    "inversion_precursor_count": prior_types.count("inversion"),
                    "duplication_precursor_count": prior_types.count("duplication"),
                    "terminal_chromosome": str(terminal_row.chromosome),
                    "terminal_pair_weight_source": str(terminal_row.pair_weight_source),
                    "terminal_pair_distance": terminal_metrics["distance"],
                    "terminal_pair_proposal_probability": terminal_metrics["proposal_probability"],
                    "terminal_missing_gate_genes": str(terminal_row.missing_gate_genes),
                }
                for event_type in ["inversion", "duplication"]:
                    effects = event_effects.get(event_type, [])
                    finite_ratios = [
                        float(item["proposal_log2_ratio_nonzero"])
                        for item in effects
                        if math.isfinite(float(item["proposal_log2_ratio_nonzero"]))
                    ]
                    trajectory_record[f"any_{event_type}_pair_activation"] = any(item["final_pair_activated"] for item in effects)
                    trajectory_record[f"any_{event_type}_final_pair_hazard_increase"] = any(item["final_pair_lethal_hazard_increased"] for item in effects)
                    trajectory_record[f"max_{event_type}_proposal_log2_ratio_nonzero"] = max(finite_ratios) if finite_ratios else np.nan
                trajectory_writer.writerow(trajectory_record)

            replay_audit.append(
                {
                    "scenario": scenario,
                    "model_key": model_key,
                    "alpha": alpha,
                    "p_event": p_event,
                    "logged_events": len(frame),
                    "gate_failing_trajectories": len(failed_ids),
                    "replayed_trajectories_with_nondeletion_precursors": scenario_replayed,
                    "nondeletion_precursor_events": scenario_precursors,
                    "sampled_total_hazard_events": len(sample_keys),
                    "replay_mismatches": scenario_mismatches,
                    "elapsed_seconds": time.perf_counter() - scenario_start,
                }
            )
            print(
                f"[{scenario_number}/{len(event_paths)}] {scenario}: failures={len(failed_ids)}, "
                f"replayed={scenario_replayed}, precursors={scenario_precursors}, "
                f"mismatches={scenario_mismatches}, elapsed={time.perf_counter() - scenario_start:.1f}s",
                flush=True,
            )
    finally:
        precursor_handle.close()
        trajectory_handle.close()

    audit_frame = pd.DataFrame(replay_audit)
    audit_frame.to_csv(args.output_dir / "cascade_replay_audit.csv", index=False)
    if int(audit_frame["replay_mismatches"].sum()) != 0:
        raise RuntimeError("Cascade replay produced topology or Gate mismatches; inspect cascade_replay_audit.csv.")

    sequence_rows = [
        {
            "scope": "scenario",
            "model_key": model,
            "alpha": alpha,
            "p_event": p_event,
            "history_class": history,
            "gate_failing_trajectories": count,
        }
        for (model, alpha, p_event, history), count in sorted(sequence_counts.items())
    ]
    sequence_rows.extend(
        {
            "scope": "overall",
            "model_key": "ALL",
            "alpha": np.nan,
            "p_event": np.nan,
            "history_class": history,
            "gate_failing_trajectories": count,
        }
        for history, count in sorted(sequence_overall.items())
    )
    sequence_frame = pd.DataFrame(sequence_rows)
    sequence_frame.to_csv(args.output_dir / "cascade_sequence_summary.csv", index=False)

    precursor_frame = pd.read_csv(precursor_path)
    for column in [
        "same_chromosome_as_terminal", "final_pair_activated", "proposal_probability_increased",
        "final_pair_lethal_hazard_increased", "total_hazard_sampled", "total_lethal_hazard_increased",
    ]:
        precursor_frame[column] = precursor_frame[column].astype(str).str.lower().eq("true")
    precursor_summary = summarize_precursors(precursor_frame)
    precursor_summary.to_csv(args.output_dir / "cascade_event_effect_summary.csv", index=False)

    trajectory_frame = pd.read_csv(trajectory_path)
    for column in [
        "any_inversion_pair_activation", "any_duplication_pair_activation",
        "any_inversion_final_pair_hazard_increase", "any_duplication_final_pair_hazard_increase",
    ]:
        trajectory_frame[column] = trajectory_frame[column].astype(str).str.lower().eq("true")
    trajectory_summary = summarize_trajectory_effects(trajectory_frame)
    trajectory_summary.to_csv(args.output_dir / "cascade_trajectory_effect_summary.csv", index=False)

    mechanism_summary = (
        precursor_frame.groupby(["precursor_event_type", "cascade_mechanism"], as_index=False)
        .agg(n_events=("trajectory_key", "size"), n_trajectories=("trajectory_key", "nunique"))
    )
    mechanism_summary["event_fraction_within_type"] = mechanism_summary["n_events"] / mechanism_summary.groupby(
        "precursor_event_type"
    )["n_events"].transform("sum")
    mechanism_summary.to_csv(args.output_dir / "cascade_mechanism_summary.csv", index=False)

    # Exploding every gene-event association into Python dictionaries can exceed
    # memory for the formal run. Reuse the chunked finalizer used for recovery.
    from finalize_rearrangement_cascade import finalize_gene_summary

    gene_summary = finalize_gene_summary(
        precursor_path,
        trajectory_path,
        args.inputs_dir,
        chunk_size=25_000,
    )
    gene_summary.to_csv(args.output_dir / "cascade_gene_summary.csv", index=False)

    top_records = [item[2] for item in sorted(top_events, key=lambda item: item[0], reverse=True)]
    pd.DataFrame(top_records).to_csv(args.output_dir / "cascade_top_events.csv", index=False)

    overall_events = precursor_summary[precursor_summary["scope"].eq("overall")]
    overall_trajectories = trajectory_summary[trajectory_summary["model_key"].eq("ALL")]
    report_lines = [
        "# Rearrangement cascade analysis",
        "",
        "## Scope",
        "",
        f"The analysis replayed {total_replayed:,} Gate-failing trajectories containing at least one prior inversion or duplication across {len(event_paths)} formal scenarios. It evaluated {total_precursors:,} non-deletion precursor events. All {total_failures:,} Gate-failing trajectories ended with deletion, as required by the copy-number-based Essentiality Gate.",
        "",
        "The primary event-level quantity is the state-conditioned probability of the exact lox pair that later produced the terminal Gate-failing deletion. For a stratified sample of precursor events, the analysis also enumerated all currently available deletion pairs and calculated the total one-step Gate-failing deletion hazard.",
        "",
        "## Sequence histories",
        "",
        markdown_table(
            sequence_frame[sequence_frame["scope"].eq("overall")].sort_values(
                "gate_failing_trajectories", ascending=False
            ),
            ["history_class", "gate_failing_trajectories"],
        ),
        "",
        "## Event-level cascade effects",
        "",
        markdown_table(
            overall_events,
            [
                "precursor_event_type", "n_precursor_events", "n_trajectories",
                "final_pair_activation_fraction", "proposal_probability_increase_fraction",
                "final_pair_lethal_hazard_increase_fraction",
                "median_proposal_log2_ratio_when_nonzero", "n_total_hazard_sampled",
                "total_lethal_hazard_increase_fraction_sample",
                "median_total_hazard_log2_ratio_nonzero_sample",
            ],
        ),
        "",
        "## Trajectory-level cascade effects",
        "",
        markdown_table(
            overall_trajectories,
            [
                "precursor_event_type", "n_trajectories",
                "trajectory_fraction_with_any_final_pair_hazard_increase",
                "trajectory_fraction_ci_low", "trajectory_fraction_ci_high",
                "trajectory_fraction_with_pair_activation",
                "median_of_trajectory_max_log2_probability_ratio",
            ],
        ),
        "",
        "## Interpretation boundary",
        "",
        "A terminal deletion is mathematically required for Gate failure under the current model. The cascade analysis therefore tests how prior inversions and duplications modify the accessibility, topology, and Gate consequence of later deletion paths. It does not establish that deletion is the sole biological cause of SCRaMbLE lethality.",
        "",
        "All-pair lethal-hazard calculations use a deterministic stratified sample of precursor states; exact terminal-pair calculations cover every logged inversion and duplication preceding Gate failure.",
    ]
    (args.output_dir / "CASCADE_ANALYSIS_REPORT.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    config = {
        "event_root": str(args.event_root),
        "inputs_dir": str(args.inputs_dir),
        "gate_name": args.gate_name,
        "scenario_count": len(event_paths),
        "logged_events": total_events,
        "gate_failing_trajectories": total_failures,
        "replayed_trajectories": total_replayed,
        "nondeletion_precursor_events": total_precursors,
        "deletion_probability": deletion_probability,
        "total_hazard_sample_per_scenario_event_type": args.total_hazard_sample_per_type,
        "seed": args.seed,
        "elapsed_seconds": time.perf_counter() - start_all,
    }
    (args.output_dir / "cascade_analysis_config.json").write_text(
        json.dumps(config, indent=2), encoding="utf-8"
    )
    print(json.dumps(config, indent=2), flush=True)


if __name__ == "__main__":
    main()
