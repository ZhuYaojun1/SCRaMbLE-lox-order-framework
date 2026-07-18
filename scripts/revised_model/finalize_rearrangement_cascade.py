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
EVENT_TYPES = ("inversion", "duplication")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Finalize cascade-analysis outputs from completed compressed replay tables."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPOSITORY_ROOT / "outputs" / "revised_model" / "cascade_analysis",
    )
    parser.add_argument(
        "--inputs-dir",
        type=Path,
        default=REPOSITORY_ROOT / "outputs" / "revised_model" / "inputs",
    )
    parser.add_argument("--chunk-size", type=int, default=25_000)
    parser.add_argument("--top-n", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20260718)
    return parser.parse_args()


def split_ids(value: Any) -> list[str]:
    if value is None or pd.isna(value):
        return []
    return [token for token in str(value).split(";") if token]


def as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.astype(str).str.lower().eq("true")


def markdown_table(frame: pd.DataFrame, columns: list[str]) -> str:
    selected = frame.loc[:, columns].copy()
    lines = ["| " + " | ".join(columns) + " |", "|" + "|".join(["---"] * len(columns)) + "|"]
    for record in selected.to_dict("records"):
        values: list[str] = []
        for column in columns:
            value = record[column]
            if isinstance(value, (float, np.floating)):
                values.append("NA" if not math.isfinite(float(value)) else f"{float(value):.4f}")
            else:
                values.append(str(value))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def load_gene_names(inputs_dir: Path) -> dict[str, str]:
    mapping = pd.read_csv(inputs_dir / "gene_lox_exact_mapping.csv", usecols=["sgd_gene_id", "gene_name"])
    mapping = mapping.dropna(subset=["sgd_gene_id"]).drop_duplicates("sgd_gene_id")
    return dict(zip(mapping["sgd_gene_id"].astype(str), mapping["gene_name"].fillna("").astype(str)))


def finalize_gene_summary(
    precursor_path: Path,
    trajectory_path: Path,
    inputs_dir: Path,
    chunk_size: int,
) -> pd.DataFrame:
    event_counts: Counter[tuple[str, str]] = Counter()
    activation_counts: Counter[tuple[str, str]] = Counter()
    hazard_counts: Counter[tuple[str, str]] = Counter()
    trajectory_counts: Counter[tuple[str, str]] = Counter()
    trajectory_activation_counts: Counter[tuple[str, str]] = Counter()
    trajectory_hazard_counts: Counter[tuple[str, str]] = Counter()

    event_columns = [
        "terminal_missing_gate_genes",
        "precursor_event_type",
        "final_pair_activated",
        "final_pair_lethal_hazard_increased",
    ]
    for chunk in pd.read_csv(precursor_path, usecols=event_columns, chunksize=chunk_size):
        chunk["final_pair_activated"] = as_bool(chunk["final_pair_activated"])
        chunk["final_pair_lethal_hazard_increased"] = as_bool(
            chunk["final_pair_lethal_hazard_increased"]
        )
        for row in chunk.itertuples(index=False):
            event_type = str(row.precursor_event_type)
            for gene in split_ids(row.terminal_missing_gate_genes):
                key = (gene, event_type)
                event_counts[key] += 1
                activation_counts[key] += int(row.final_pair_activated)
                hazard_counts[key] += int(row.final_pair_lethal_hazard_increased)

    trajectory_columns = [
        "terminal_missing_gate_genes",
        "inversion_precursor_count",
        "duplication_precursor_count",
        "any_inversion_pair_activation",
        "any_duplication_pair_activation",
        "any_inversion_final_pair_hazard_increase",
        "any_duplication_final_pair_hazard_increase",
    ]
    for chunk in pd.read_csv(trajectory_path, usecols=trajectory_columns, chunksize=chunk_size):
        for event_type in EVENT_TYPES:
            chunk[f"any_{event_type}_pair_activation"] = as_bool(
                chunk[f"any_{event_type}_pair_activation"]
            )
            chunk[f"any_{event_type}_final_pair_hazard_increase"] = as_bool(
                chunk[f"any_{event_type}_final_pair_hazard_increase"]
            )
        for row in chunk.itertuples(index=False):
            genes = split_ids(row.terminal_missing_gate_genes)
            for event_type in EVENT_TYPES:
                if int(getattr(row, f"{event_type}_precursor_count")) <= 0:
                    continue
                activated = bool(getattr(row, f"any_{event_type}_pair_activation"))
                increased = bool(getattr(row, f"any_{event_type}_final_pair_hazard_increase"))
                for gene in genes:
                    key = (gene, event_type)
                    trajectory_counts[key] += 1
                    trajectory_activation_counts[key] += int(activated)
                    trajectory_hazard_counts[key] += int(increased)

    gene_names = load_gene_names(inputs_dir)
    rows: list[dict[str, Any]] = []
    for gene, event_type in sorted(event_counts):
        n_events = event_counts[(gene, event_type)]
        n_trajectories = trajectory_counts[(gene, event_type)]
        rows.append(
            {
                "sgd_gene_id": gene,
                "gene_name": gene_names.get(gene, ""),
                "precursor_event_type": event_type,
                "n_precursor_events": n_events,
                "n_trajectories": n_trajectories,
                "pair_activation_events": activation_counts[(gene, event_type)],
                "hazard_increase_events": hazard_counts[(gene, event_type)],
                "hazard_increase_event_fraction": hazard_counts[(gene, event_type)] / n_events,
                "trajectories_with_pair_activation": trajectory_activation_counts[(gene, event_type)],
                "trajectories_with_hazard_increase": trajectory_hazard_counts[(gene, event_type)],
                "hazard_increase_trajectory_fraction": (
                    trajectory_hazard_counts[(gene, event_type)] / n_trajectories
                    if n_trajectories
                    else np.nan
                ),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["hazard_increase_events", "trajectories_with_hazard_increase", "n_trajectories"],
        ascending=False,
    )


def finalize_event_derivatives(
    precursor_path: Path,
    chunk_size: int,
    top_n: int,
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    top_chunks: list[pd.DataFrame] = []
    amplification_chunks: list[pd.DataFrame] = []
    activation_chunks: list[pd.DataFrame] = []
    lethality_activation_chunks: list[pd.DataFrame] = []
    total_hazard_increase_chunks: list[pd.DataFrame] = []
    sampled_chunks: list[pd.DataFrame] = []
    chromosome_chunks: list[pd.DataFrame] = []
    timing_chunks: list[pd.DataFrame] = []
    amplification_counts: Counter[tuple[str, str]] = Counter()
    amplification_trajectories: dict[tuple[str, str], set[str]] = {}
    thresholds = {
        "at_least_2_fold": 1.0,
        "at_least_10_fold": math.log2(10.0),
        "at_least_100_fold": math.log2(100.0),
    }

    for chunk in pd.read_csv(precursor_path, chunksize=chunk_size):
        for column in [
            "same_chromosome_as_terminal",
            "final_pair_activated",
            "final_pair_lethal_hazard_increased",
            "total_hazard_sampled",
            "total_lethal_hazard_increased",
        ]:
            chunk[column] = as_bool(chunk[column])

        ratio = pd.to_numeric(chunk["proposal_log2_ratio_nonzero"], errors="coerce")
        probability_after = pd.to_numeric(chunk["proposal_probability_after"], errors="coerce")
        chunk["cascade_score"] = ratio
        activated = chunk["final_pair_activated"]
        chunk.loc[activated, "cascade_score"] = 1_000.0 + np.log10(
            probability_after.loc[activated].clip(lower=1e-300)
        )
        finite = chunk[np.isfinite(chunk["cascade_score"])].nlargest(top_n, "cascade_score")
        top_chunks.append(finite)
        amplification_chunks.append(
            chunk[
                (~chunk["final_pair_activated"])
                & pd.to_numeric(chunk["proposal_log2_ratio_nonzero"], errors="coerce").gt(0)
            ].nlargest(top_n, "proposal_log2_ratio_nonzero")
        )
        activation_chunks.append(
            chunk[chunk["final_pair_activated"]].nlargest(top_n, "proposal_probability_after")
        )
        lethality_activation_chunks.append(
            chunk[as_bool(chunk["final_pair_lethality_activated"])].nlargest(
                top_n, "proposal_probability_after"
            )
        )
        for event_type in EVENT_TYPES:
            event_type_mask = chunk["precursor_event_type"].eq(event_type)
            for threshold_name, threshold_value in thresholds.items():
                selected = chunk[event_type_mask & ratio.ge(threshold_value)]
                key = (event_type, threshold_name)
                amplification_counts[key] += len(selected)
                amplification_trajectories.setdefault(key, set()).update(
                    selected["trajectory_key"].astype(str)
                )
            activated_rows = chunk[event_type_mask & chunk["final_pair_activated"]]
            activation_key = (event_type, "new_pair_activation")
            amplification_counts[activation_key] += len(activated_rows)
            amplification_trajectories.setdefault(activation_key, set()).update(
                activated_rows["trajectory_key"].astype(str)
            )

        sampled = chunk[chunk["total_hazard_sampled"]].copy()
        if len(sampled):
            sampled_chunks.append(sampled)
            total_hazard_increase_chunks.append(
                sampled[
                    pd.to_numeric(sampled["total_hazard_log2_ratio_nonzero"], errors="coerce").gt(0)
                ].nlargest(top_n, "total_hazard_log2_ratio_nonzero")
            )

        chromosome_chunks.append(
            chunk.groupby(
                ["precursor_event_type", "precursor_chromosome", "terminal_chromosome"],
                as_index=False,
            ).agg(
                n_precursor_events=("trajectory_key", "size"),
                pair_activation_events=("final_pair_activated", "sum"),
                final_pair_hazard_increase_events=("final_pair_lethal_hazard_increased", "sum"),
            )
        )

        timing = chunk.assign(
            timing_class=np.select(
                [chunk["events_until_failure"].eq(1), chunk["events_until_failure"].eq(2)],
                ["immediate_predecessor", "two_events_before_failure"],
                default="three_or_more_events_before_failure",
            )
        )
        timing_chunks.append(
            timing.groupby(["precursor_event_type", "timing_class"], as_index=False).agg(
                n_precursor_events=("trajectory_key", "size"),
                final_pair_hazard_increase_events=("final_pair_lethal_hazard_increased", "sum"),
            )
        )

    top_events = pd.concat(top_chunks, ignore_index=True).nlargest(top_n, "cascade_score")
    event_lists = {
        "cascade_top_events.csv": top_events,
        "cascade_top_probability_amplifications.csv": pd.concat(
            amplification_chunks, ignore_index=True
        ).nlargest(top_n, "proposal_log2_ratio_nonzero"),
        "cascade_top_pair_activations.csv": pd.concat(
            activation_chunks, ignore_index=True
        ).nlargest(top_n, "proposal_probability_after"),
        "cascade_top_lethality_activations.csv": pd.concat(
            lethality_activation_chunks, ignore_index=True
        ).nlargest(top_n, "proposal_probability_after"),
        "cascade_top_total_hazard_increases.csv": pd.concat(
            total_hazard_increase_chunks, ignore_index=True
        ).nlargest(top_n, "total_hazard_log2_ratio_nonzero"),
    }
    for frame in event_lists.values():
        frame["proposal_probability_fold_change_nonzero"] = np.power(
            2.0, pd.to_numeric(frame["proposal_log2_ratio_nonzero"], errors="coerce")
        )

    sampled_frame = pd.concat(sampled_chunks, ignore_index=True)
    sample_groups: list[tuple[str, pd.DataFrame]] = [("ALL", sampled_frame)]
    sample_groups.extend(list(sampled_frame.groupby("model_key", sort=True)))
    sample_rows: list[dict[str, Any]] = []
    for model_key, model_group in sample_groups:
        for event_type, group in model_group.groupby("precursor_event_type", sort=True):
            ratio = pd.to_numeric(group["total_hazard_log2_ratio_nonzero"], errors="coerce").dropna()
            before = pd.to_numeric(group["total_lethal_hazard_before"], errors="coerce")
            after = pd.to_numeric(group["total_lethal_hazard_after"], errors="coerce")
            delta = after - before
            tolerance = 1e-18
            increased = as_bool(group["total_lethal_hazard_increased"])
            decreased = delta.lt(-tolerance)
            unchanged = ~(increased | decreased)
            sample_rows.append(
                {
                    "model_key": model_key,
                    "precursor_event_type": event_type,
                    "n_sampled_events": len(group),
                    "hazard_increase_events": int(increased.sum()),
                    "hazard_decrease_events": int(decreased.sum()),
                    "hazard_unchanged_events": int(unchanged.sum()),
                    "hazard_increase_fraction": float(increased.mean()),
                    "hazard_decrease_fraction": float(decreased.mean()),
                    "median_hazard_log2_ratio_nonzero": float(ratio.median()) if len(ratio) else np.nan,
                    "q025_hazard_log2_ratio_nonzero": float(ratio.quantile(0.025)) if len(ratio) else np.nan,
                    "q975_hazard_log2_ratio_nonzero": float(ratio.quantile(0.975)) if len(ratio) else np.nan,
                }
            )
    total_hazard_summary = pd.DataFrame(sample_rows)

    chromosome_summary = (
        pd.concat(chromosome_chunks, ignore_index=True)
        .groupby(["precursor_event_type", "precursor_chromosome", "terminal_chromosome"], as_index=False)
        .sum(numeric_only=True)
    )
    chromosome_summary["pair_activation_fraction"] = (
        chromosome_summary["pair_activation_events"] / chromosome_summary["n_precursor_events"]
    )
    chromosome_summary["final_pair_hazard_increase_fraction"] = (
        chromosome_summary["final_pair_hazard_increase_events"]
        / chromosome_summary["n_precursor_events"]
    )

    timing_summary = (
        pd.concat(timing_chunks, ignore_index=True)
        .groupby(["precursor_event_type", "timing_class"], as_index=False)
        .sum(numeric_only=True)
    )
    timing_summary["final_pair_hazard_increase_fraction"] = (
        timing_summary["final_pair_hazard_increase_events"] / timing_summary["n_precursor_events"]
    )
    amplification_rows = []
    for event_type in EVENT_TYPES:
        for threshold_name in [*thresholds, "new_pair_activation"]:
            key = (event_type, threshold_name)
            amplification_rows.append(
                {
                    "precursor_event_type": event_type,
                    "change_class": threshold_name,
                    "n_precursor_events": amplification_counts[key],
                    "n_trajectories": len(amplification_trajectories.get(key, set())),
                }
            )
    amplification_summary = pd.DataFrame(amplification_rows)
    return event_lists, total_hazard_summary, chromosome_summary, timing_summary, amplification_summary


def write_report(
    output_dir: Path,
    audit: pd.DataFrame,
    sequence: pd.DataFrame,
    event_summary: pd.DataFrame,
    trajectory_summary: pd.DataFrame,
    mechanism_summary: pd.DataFrame,
    total_hazard_summary: pd.DataFrame,
    amplification_summary: pd.DataFrame,
    gene_summary: pd.DataFrame,
) -> None:
    overall_events = event_summary[event_summary["scope"].eq("overall")].copy()
    overall_trajectories = trajectory_summary[trajectory_summary["model_key"].eq("ALL")].copy()
    overall_sequence = sequence[sequence["scope"].eq("overall")].copy()
    overall_hazard = total_hazard_summary[total_hazard_summary["model_key"].eq("ALL")].copy()
    total_failures = int(audit["gate_failing_trajectories"].sum())
    replayed = int(audit["replayed_trajectories_with_nondeletion_precursors"].sum())
    precursors = int(audit["nondeletion_precursor_events"].sum())
    failure_fraction = replayed / total_failures
    top_genes = gene_summary.head(15).copy()

    report = [
        "# Rearrangement cascade analysis",
        "",
        "## Principal finding",
        "",
        f"Among {total_failures:,} Essentiality-Gate-failing trajectories, {replayed:,} ({failure_fraction:.2%}) contained at least one inversion or duplication before the terminal deletion. Exact state replay evaluated {precursors:,} such precursor events with no topology, copy-number, or Gate-status mismatches. Prior rearrangements therefore occurred in a substantial fraction of terminal deletion histories, although their effects were heterogeneous rather than uniformly risk enhancing.",
        "",
        "The copy-number Gate makes deletion the only event capable of producing an immediate Gate failure. The analysis consequently asks a narrower and testable question: whether earlier inversion or duplication changed the accessibility, proposal probability, or Gate consequence of the deletion path that later ended the trajectory.",
        "",
        "## Failure-history composition",
        "",
        markdown_table(overall_sequence.sort_values("gate_failing_trajectories", ascending=False), ["history_class", "gate_failing_trajectories"]),
        "",
        "## Exact terminal-pair effects",
        "",
        markdown_table(
            overall_events,
            [
                "precursor_event_type",
                "n_precursor_events",
                "n_trajectories",
                "final_pair_activation_fraction",
                "proposal_probability_increase_fraction",
                "final_pair_lethal_hazard_increase_fraction",
                "median_proposal_log2_ratio_when_nonzero",
            ],
        ),
        "",
        markdown_table(
            overall_trajectories,
            [
                "precursor_event_type",
                "n_trajectories",
                "trajectory_fraction_with_any_final_pair_hazard_increase",
                "trajectory_fraction_ci_low",
                "trajectory_fraction_ci_high",
                "trajectory_fraction_with_pair_activation",
            ],
        ),
        "",
        "Inversions most often left the exact future terminal pair unchanged, but a minority shortened its current lox-order distance, changed its sampling weight, or converted its interval into a Gate-failing deletion span. Duplications could create a new copy-specific terminal pair, but more often expanded the competing pair space or buffered essential-gene copy number.",
        "",
        "## Large terminal-pair probability changes",
        "",
        markdown_table(
            amplification_summary,
            ["precursor_event_type", "change_class", "n_precursor_events", "n_trajectories"],
        ),
        "",
        "Finite fold-change classes exclude newly activated pairs whose pre-event probability was zero. New-pair activation is therefore reported separately.",
        "",
        "## Prospective all-pair Gate-failing deletion hazard",
        "",
        markdown_table(
            overall_hazard,
            [
                "precursor_event_type",
                "n_sampled_events",
                "hazard_increase_events",
                "hazard_decrease_events",
                "hazard_unchanged_events",
                "hazard_increase_fraction",
                "hazard_decrease_fraction",
                "median_hazard_log2_ratio_nonzero",
            ],
        ),
        "",
        "These all-pair values come from a deterministic, scenario-stratified sample of up to 50 precursor events of each type per scenario. They enumerate every active deletion pair in the state immediately before and after the precursor event and therefore summarize the one-step Gate-failing deletion hazard, not only the pair eventually observed in that same trajectory.",
        "",
        "## Mechanistic categories",
        "",
        markdown_table(
            mechanism_summary.sort_values(["precursor_event_type", "n_events"], ascending=[True, False]),
            ["precursor_event_type", "cascade_mechanism", "n_events", "n_trajectories", "event_fraction_within_type"],
        ),
        "",
        "## Missing-gene associations",
        "",
        markdown_table(
            top_genes,
            [
                "sgd_gene_id",
                "gene_name",
                "precursor_event_type",
                "n_trajectories",
                "hazard_increase_events",
                "hazard_increase_event_fraction",
            ],
        ),
        "",
        "Gene rows describe associations with the essential genes absent after the later terminal deletion. They do not show that the precursor event deleted those genes or identify a biological lethal mechanism outside the modeled Essentiality Gate.",
        "",
        "## Interpretation boundary",
        "",
        "The exact terminal-pair analysis is conditioned on trajectories that later failed the Gate and is therefore retrospective; it quantifies path modification but does not estimate a population causal effect. The sampled all-pair analysis provides a broader state-conditioned risk measure, but it remains a one-step model hazard under the fitted event and pair-sampling rules. The results support a rearrangement-cascade interpretation in which earlier events can reshape later deletion accessibility, while also showing that most individual precursor events do not increase the eventual terminal-pair hazard.",
    ]
    (output_dir / "CASCADE_ANALYSIS_REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def validate_outputs(
    output_dir: Path,
    precursor_path: Path,
    trajectory_path: Path,
    audit: pd.DataFrame,
    sequence: pd.DataFrame,
    event_summary: pd.DataFrame,
    mechanism_summary: pd.DataFrame,
    total_hazard_summary: pd.DataFrame,
    chunk_size: int,
) -> None:
    precursor_rows = 0
    probability_bound_violations = 0
    activation_violations = 0
    lethality_activation_violations = 0
    scenario_precursors: Counter[str] = Counter()
    validation_columns = [
        "scenario",
        "proposal_probability_before",
        "proposal_probability_after",
        "final_pair_activated",
        "final_pair_exists_before",
        "final_pair_exists_after",
        "final_pair_lethality_activated",
        "final_pair_lethal_before",
        "final_pair_lethal_after",
    ]
    for chunk in pd.read_csv(precursor_path, usecols=validation_columns, chunksize=chunk_size):
        precursor_rows += len(chunk)
        scenario_precursors.update(chunk["scenario"].astype(str).value_counts().to_dict())
        before = pd.to_numeric(chunk["proposal_probability_before"], errors="coerce")
        after = pd.to_numeric(chunk["proposal_probability_after"], errors="coerce")
        probability_bound_violations += int((~before.between(0.0, 1.0)).sum())
        probability_bound_violations += int((~after.between(0.0, 1.0)).sum())
        activated = as_bool(chunk["final_pair_activated"])
        exists_before = as_bool(chunk["final_pair_exists_before"])
        exists_after = as_bool(chunk["final_pair_exists_after"])
        activation_violations += int((activated & (exists_before | ~exists_after)).sum())
        lethality_activated = as_bool(chunk["final_pair_lethality_activated"])
        lethal_before = as_bool(chunk["final_pair_lethal_before"])
        lethal_after = as_bool(chunk["final_pair_lethal_after"])
        lethality_activation_violations += int(
            (lethality_activated & (lethal_before | ~lethal_after)).sum()
        )

    trajectory_rows = 0
    scenario_trajectories: Counter[str] = Counter()
    for chunk in pd.read_csv(
        trajectory_path, usecols=["scenario", "trajectory_key"], chunksize=chunk_size
    ):
        trajectory_rows += len(chunk)
        scenario_trajectories.update(chunk["scenario"].astype(str).value_counts().to_dict())

    audit_precursors = dict(
        zip(audit["scenario"].astype(str), audit["nondeletion_precursor_events"].astype(int))
    )
    audit_trajectories = dict(
        zip(
            audit["scenario"].astype(str),
            audit["replayed_trajectories_with_nondeletion_precursors"].astype(int),
        )
    )
    overall_event_counts = (
        event_summary[event_summary["scope"].eq("overall")]
        .set_index("precursor_event_type")["n_precursor_events"]
        .astype(int)
        .to_dict()
    )
    mechanism_counts = (
        mechanism_summary.groupby("precursor_event_type")["n_events"].sum().astype(int).to_dict()
    )
    overall_hazard = total_hazard_summary[total_hazard_summary["model_key"].eq("ALL")]
    hazard_partition_ok = bool(
        (
            overall_hazard["hazard_increase_events"]
            + overall_hazard["hazard_decrease_events"]
            + overall_hazard["hazard_unchanged_events"]
            == overall_hazard["n_sampled_events"]
        ).all()
    )
    checks = [
        ("All 88 scenarios were audited", len(audit) == 88, len(audit)),
        ("Replay mismatches are zero", int(audit["replay_mismatches"].sum()) == 0, int(audit["replay_mismatches"].sum())),
        ("Precursor detail row count matches audit", precursor_rows == int(audit["nondeletion_precursor_events"].sum()), precursor_rows),
        ("Precursor counts match audit scenario by scenario", dict(scenario_precursors) == audit_precursors, "matched" if dict(scenario_precursors) == audit_precursors else "mismatch"),
        ("Trajectory detail row count matches audit", trajectory_rows == int(audit["replayed_trajectories_with_nondeletion_precursors"].sum()), trajectory_rows),
        ("Trajectory counts match audit scenario by scenario", dict(scenario_trajectories) == audit_trajectories, "matched" if dict(scenario_trajectories) == audit_trajectories else "mismatch"),
        ("Failure-history categories close to all Gate failures", int(sequence[sequence["scope"].eq("overall")]["gate_failing_trajectories"].sum()) == int(audit["gate_failing_trajectories"].sum()), int(sequence[sequence["scope"].eq("overall")]["gate_failing_trajectories"].sum())),
        ("Event summary and mechanism categories agree", overall_event_counts == mechanism_counts, str(overall_event_counts)),
        ("Total-hazard sample count matches audit", int(overall_hazard["n_sampled_events"].sum()) == int(audit["sampled_total_hazard_events"].sum()), int(overall_hazard["n_sampled_events"].sum())),
        ("Total-hazard direction categories are exhaustive", hazard_partition_ok, "matched" if hazard_partition_ok else "mismatch"),
        ("All proposal probabilities are within [0, 1]", probability_bound_violations == 0, probability_bound_violations),
        ("Pair-activation flags agree with pair existence", activation_violations == 0, activation_violations),
        ("Lethality-activation flags agree with Gate status", lethality_activation_violations == 0, lethality_activation_violations),
    ]
    failed = [name for name, passed, _ in checks if not passed]
    lines = [
        "# Cascade analysis validation",
        "",
        f"Overall status: {'PASS' if not failed else 'FAIL'}",
        "",
        "| Check | Status | Observed |",
        "|---|---|---|",
    ]
    lines.extend(
        f"| {name} | {'PASS' if passed else 'FAIL'} | {observed} |"
        for name, passed, observed in checks
    )
    lines.extend(
        [
            "",
            "The validation checks accounting closure, scenario-level agreement, categorical exhaustiveness, probability bounds, and logical consistency of activation flags. It does not convert the retrospective terminal-pair analysis into a causal population estimate.",
        ]
    )
    (output_dir / "CASCADE_ANALYSIS_VALIDATION.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    if failed:
        raise RuntimeError("Cascade validation failed: " + "; ".join(failed))


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    precursor_path = output_dir / "cascade_precursor_events.csv.xz"
    trajectory_path = output_dir / "cascade_trajectory_effects.csv.xz"
    required = [
        precursor_path,
        trajectory_path,
        output_dir / "cascade_replay_audit.csv",
        output_dir / "cascade_sequence_summary.csv",
        output_dir / "cascade_event_effect_summary.csv",
        output_dir / "cascade_trajectory_effect_summary.csv",
        output_dir / "cascade_mechanism_summary.csv",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing completed replay outputs: " + "; ".join(missing))

    print("Streaming gene-level aggregation...", flush=True)
    gene_summary = finalize_gene_summary(
        precursor_path, trajectory_path, args.inputs_dir.resolve(), args.chunk_size
    )
    gene_summary.to_csv(output_dir / "cascade_gene_summary.csv", index=False)

    print("Streaming top-event, chromosome, timing, and total-hazard summaries...", flush=True)
    event_lists, total_hazard, chromosome, timing, amplification = finalize_event_derivatives(
        precursor_path, args.chunk_size, args.top_n
    )
    for filename, frame in event_lists.items():
        frame.to_csv(output_dir / filename, index=False)
    total_hazard.to_csv(output_dir / "cascade_total_hazard_sample_summary.csv", index=False)
    chromosome.to_csv(output_dir / "cascade_chromosome_summary.csv", index=False)
    timing.to_csv(output_dir / "cascade_timing_summary.csv", index=False)
    amplification.to_csv(output_dir / "cascade_amplification_threshold_summary.csv", index=False)

    audit = pd.read_csv(output_dir / "cascade_replay_audit.csv")
    sequence = pd.read_csv(output_dir / "cascade_sequence_summary.csv")
    event_summary = pd.read_csv(output_dir / "cascade_event_effect_summary.csv")
    trajectory_summary = pd.read_csv(output_dir / "cascade_trajectory_effect_summary.csv")
    mechanism_summary = pd.read_csv(output_dir / "cascade_mechanism_summary.csv")
    write_report(
        output_dir,
        audit,
        sequence,
        event_summary,
        trajectory_summary,
        mechanism_summary,
        total_hazard,
        amplification,
        gene_summary,
    )
    validate_outputs(
        output_dir,
        precursor_path,
        trajectory_path,
        audit,
        sequence,
        event_summary,
        mechanism_summary,
        total_hazard,
        args.chunk_size,
    )

    config = {
        "source": "completed deterministic replay outputs",
        "scenario_count": int(len(audit)),
        "logged_events": int(audit["logged_events"].sum()),
        "gate_failing_trajectories": int(audit["gate_failing_trajectories"].sum()),
        "replayed_trajectories_with_nondeletion_precursors": int(
            audit["replayed_trajectories_with_nondeletion_precursors"].sum()
        ),
        "nondeletion_precursor_events": int(audit["nondeletion_precursor_events"].sum()),
        "sampled_total_hazard_events": int(audit["sampled_total_hazard_events"].sum()),
        "replay_mismatches": int(audit["replay_mismatches"].sum()),
        "postprocessing_chunk_size": args.chunk_size,
        "top_event_count": args.top_n,
        "selection_seed": args.seed,
    }
    (output_dir / "cascade_analysis_config.json").write_text(
        json.dumps(config, indent=2), encoding="utf-8"
    )
    print(json.dumps(config, indent=2), flush=True)


if __name__ == "__main__":
    main()
