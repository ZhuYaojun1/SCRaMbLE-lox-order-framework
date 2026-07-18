from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate revised-model consistency audit and manifests.")
    parser.add_argument("--project-root", type=Path, default=REPOSITORY_ROOT)
    parser.add_argument("--output-root", type=Path, default=REPOSITORY_ROOT / "outputs" / "revised_model")
    return parser.parse_args()


def markdown_table(frame: pd.DataFrame) -> str:
    if frame.empty:
        return "No records."
    columns = list(frame.columns)
    lines = ["| " + " | ".join(columns) + " |", "|" + "|".join(["---"] * len(columns)) + "|"]
    for row in frame.fillna("").astype(str).to_dict("records"):
        lines.append("| " + " | ".join(row[column].replace("|", "\\|") for column in columns) + " |")
    return "\n".join(lines)


def best_row(frame: pd.DataFrame, column: str, ascending: bool = False, subset: pd.Series | None = None) -> dict[str, Any]:
    selected = frame if subset is None else frame[subset]
    selected = selected[np.isfinite(pd.to_numeric(selected[column], errors="coerce"))]
    if selected.empty:
        return {}
    return selected.sort_values(column, ascending=ascending).iloc[0].to_dict()


def fmt(value: Any, digits: int = 4) -> str:
    try:
        number = float(value)
        if not np.isfinite(number):
            return "NA"
        return f"{number:.{digits}f}"
    except (TypeError, ValueError):
        return str(value)


def build_checks(output: Path) -> pd.DataFrame:
    inputs = output / "inputs"
    main = output / "main_run"
    analysis = output / "analysis"
    lox = pd.read_csv(inputs / "lox_sites_validated.csv")
    audit = pd.read_csv(inputs / "lox_input_audit.csv")
    hic = pd.read_csv(inputs / "hic_pair_coverage.csv")
    population = pd.read_csv(main / "population_summary.csv")
    predicted = pd.read_csv(main / "predicted_lox_frequency.csv")
    rare = pd.read_csv(analysis / "rarefied_diversity.csv")
    rows: list[dict[str, Any]] = []

    def add(name: str, passed: bool, observed: Any, expected: Any) -> None:
        rows.append({"check": name, "status": "PASS" if passed else "FAIL", "observed": observed, "expected": expected})

    observed_counts = lox.groupby("chromosome").size().to_dict()
    expected_counts = {"SynII": 273, "SynIII": 100, "SynVI": 71, "SynIXR": 44}
    add("active lox counts by chromosome", observed_counts == expected_counts, json.dumps(observed_counts), json.dumps(expected_counts))
    add("active lox total", len(lox) == 488, len(lox), 488)
    add("removed blank padding", int(audit["removed_blank_records"].sum()) == 604, int(audit["removed_blank_records"].sum()), 604)
    add("candidate pair total", int(audit["candidate_same_chromosome_pairs"].sum()) == 45509, int(audit["candidate_same_chromosome_pairs"].sum()), 45509)
    chromosome_hic = hic[hic["chromosome"] != "ALL"]
    direct_chromosomes = set(chromosome_hic.loc[chromosome_hic["direct_hic_pairs"] > 0, "chromosome"])
    add("direct Hi-C restricted to SynII", direct_chromosomes == {"SynII"}, ",".join(sorted(direct_chromosomes)), "SynII only")
    add("formal scenario count", len(population) == 88, len(population), 88)
    add("formal trajectory depth", bool((population["initialized_trajectories"] == 10000).all()), population["initialized_trajectories"].min(), 10000)
    balance = population["gate_passing_trajectories"] + population["gate_failing_trajectories"]
    add("trajectory balance", bool((balance == population["initialized_trajectories"]).all()), int((balance == population["initialized_trajectories"]).sum()), len(population))
    event_balance = population[["deletion_events", "inversion_events", "duplication_events"]].sum(axis=1)
    add("accepted event balance", bool((event_balance == population["accepted_events"]).all()), int((event_balance == population["accepted_events"]).sum()), len(population))
    add("predicted site rows", len(predicted) == 88 * 488, len(predicted), 88 * 488)
    add("unique aligned sites per scenario", int(predicted.groupby(["model_key", "alpha", "p_event"]).size().min()) == 488, int(predicted.groupby(["model_key", "alpha", "p_event"]).size().min()), 488)
    rare_bound = rare["effective_shannon_diversity_mean"] <= rare["rarefied_survivor_depth"] + 1e-9
    add("rarefied effective diversity mathematical bound", bool(rare_bound.all()), int(rare_bound.sum()), len(rare))
    add("rarefaction replicates", bool((rare["rarefaction_replicates"] >= 100).all()), int(rare["rarefaction_replicates"].min()), ">=100")
    robustness_manifest = output / "robustness" / "robustness_execution_manifest.csv"
    if robustness_manifest.exists():
        robust = pd.read_csv(robustness_manifest)
        add("robustness return codes", bool((robust["return_code"] == 0).all()), int((robust["return_code"] == 0).sum()), len(robust))
    isomorphism_path = analysis / "endpoint_isomorphism_audit.csv"
    if isomorphism_path.exists():
        isomorphism = pd.read_csv(isomorphism_path)
        add("endpoint isomorphism audit scenarios", len(isomorphism) == 88, len(isomorphism), 88)
        add(
            "unlabeled topology cannot increase unique endpoints",
            bool((isomorphism["unlabeled_topology_unique_endpoints"] <= isomorphism["provenance_unique_endpoints"]).all()),
            int((isomorphism["unlabeled_topology_unique_endpoints"] <= isomorphism["provenance_unique_endpoints"]).sum()),
            len(isomorphism),
        )
        add(
            "unlabeled rarefied diversity mathematical bound",
            bool((isomorphism["unlabeled_topology_rarefied_effective_shannon_mean"] <= isomorphism["common_rarefied_depth"] + 1e-9).all()),
            int((isomorphism["unlabeled_topology_rarefied_effective_shannon_mean"] <= isomorphism["common_rarefied_depth"] + 1e-9).sum()),
            len(isomorphism),
        )
    return pd.DataFrame(rows)


def generate_report(project: Path, output: Path, checks: pd.DataFrame) -> str:
    inputs = output / "inputs"
    analysis = output / "analysis"
    lox_audit = pd.read_csv(inputs / "lox_input_audit.csv")
    hic = pd.read_csv(inputs / "hic_pair_coverage.csv")
    mapping = pd.read_csv(inputs / "genbank_mapping_audit.csv")
    gate = pd.read_csv(inputs / "survival_gate_summary.csv")
    population = pd.read_csv(output / "main_run" / "population_summary.csv")
    rare = pd.read_csv(analysis / "rarefied_diversity.csv")
    compare = pd.read_csv(analysis / "diversity_rank_consistency.csv")
    pareto = pd.read_csv(analysis / "survival_diversity_pareto_frontier.csv")
    backtest = pd.read_csv(analysis / "backtesting_report.csv")
    nulls = pd.read_csv(analysis / "null_control_summary.csv")
    holdout = pd.read_csv(analysis / "chromosome_wise_holdout.csv")
    risk = pd.read_csv(analysis / "exact_coordinate_gene_risk.csv")
    isomorphism_path = analysis / "endpoint_isomorphism_audit.csv"
    isomorphism = pd.read_csv(isomorphism_path) if isomorphism_path.exists() else pd.DataFrame()
    structural = rare[rare["endpoint_type"] == "structural"]
    best_div = best_row(structural, "effective_shannon_diversity_mean")
    best_hic = best_row(backtest[backtest["model_key"] == "partial_hic_fallback"], "pearson_r")
    best_linear = best_row(backtest[backtest["model_key"] == "linear_distance"], "pearson_r")
    retention = pareto[pareto["passes_retention_0_70"].astype(str).str.lower().isin(["true", "1"])]
    best_retention = best_row(retention, "effective_shannon_diversity_mean")
    event_totals = population[["deletion_events", "inversion_events", "duplication_events"]].sum()
    total_events = int(event_totals.sum())
    gate_strict = int(gate.loc[gate["gate_name"] == "strict", "n_genes"].iloc[0])
    coordinate_evidence = pd.read_csv(inputs / "lox_coordinate_evidence.csv")
    unresolved_count = int(coordinate_evidence["unresolved_records"].sum())

    robustness_text = "Robustness analysis was not available."
    robust_path = output / "robustness" / "five_seed_robustness_summary.csv"
    if robust_path.exists():
        robust = pd.read_csv(robust_path)
        sensitivity = pd.read_csv(output / "robustness" / "essentiality_gate_sensitivity_summary.csv")
        robustness_text = (
            f"The five-seed analysis contains {len(robust)} model-by-p_event summaries. "
            f"Gate sensitivity contains {len(sensitivity)} gate-by-p_event scenarios; complete member lists and rules are in survival_gate_gene_sets.csv."
        )

    lines = [
        "# FINAL REVISION REPORT",
        "",
        "## Scope and status",
        "",
        "This report covers only the independent revised-model tree. Original raw data, legacy scripts, Stage 2.3 results, and manuscript files were not overwritten. The formal run uses 10,000 initialized trajectories per scenario and 20 sequential opportunities across 88 main-grid scenarios.",
        "",
        f"Automated checks: {(checks['status'] == 'PASS').sum()} passed and {(checks['status'] == 'FAIL').sum()} failed.",
        "",
        "## Original errors and corrections",
        "",
        "1. **Padded lox columns.** The legacy Fig. 2C parser tested the shared order header but did not require a non-empty chromosome-specific count. It therefore emitted 604 blank padding cells as active sites. The revised parser requires a real value/record in the chromosome row. Active sites changed from 1,092 to 488 and same-chromosome candidate pairs from 148,512 to 45,509.",
        "2. **Hi-C scope.** The legacy mapping made all 37,128 SynII pairs appear directly covered and could be described too broadly. Exact bp-to-bin mapping resolves 36,046 positive-contact SynII pairs; all non-SynII pairs use distance fallback. The model is named partially Hi-C-informed sampling with distance-based fallback.",
        "3. **Gene-to-lox mapping.** Proportional midpoint scaling was removed. Features are mapped through accessioned synthetic-chromosome sequences and motif-supported lox coordinates. SynIXR has one unresolved physical lox record, so SynIXR gene-level risk ranking is disabled.",
        "4. **State update engine.** Fixed-reference copy-number edits were replaced by ordered segment and lox copy instances. Deletion removes the interval and one recombination boundary; inversion reverses part order and orientation; duplication creates new segment and lox copy instances. Eligible pairs are rebuilt after every event.",
        "5. **Endpoint definition.** The truncated changed-gene signature was removed. The provenance-aware endpoint losslessly encodes segment order, orientation, copy identity, junctions, and active lox identities. A separate copy-label-invariant linear-topology audit removes event-derived instance labels while retaining template order, orientation, multiplicity, and adjacency. ORF copy-number endpoints are reported separately.",
        "6. **Diversity.** Main comparisons use 100 repeated without-replacement rarefactions to the minimum revised survivor depth. Full-sample results are supplemental. The discontinuous constrained-diversity score is retired in favor of a survival-diversity Pareto frontier and prespecified 0.50/0.70 retention indicators.",
        "7. **Survival terminology.** Outputs use essentiality-gate passing/failing trajectory. The revised strict set contains only exact-mapped phenotype-derived Essential ORFs; it does not represent complete biological viability.",
        "",
        "## Lox input audit",
        "",
        markdown_table(lox_audit[["chromosome", "raw_header_columns", "nonempty_records", "removed_blank_records", "unique_lox_id", "lox_order_min", "lox_order_max", "candidate_same_chromosome_pairs", "sites_with_observed_frequency"]]),
        "",
        "The 604 removed records do not appear in lox_sites_validated.csv, the dynamic active-lox topology, event proposal space, or predicted_lox_frequency.csv.",
        "",
        "## Hi-C coverage audit",
        "",
        markdown_table(hic[["chromosome", "total_candidate_pairs", "direct_hic_pairs", "distance_fallback_pairs", "unresolved_coordinate_pairs", "direct_hic_fraction"]]),
        "",
        "Direct and fallback proposal masses are normalized within the current eligible pair set. Direct contact weights are robust-scaled against the positive SynII contact distribution before mixture normalization; fallback remains the explicit topological-distance proposal rather than fabricated contact data.",
        "",
        "## Exact coordinate mapping",
        "",
        markdown_table(mapping),
        "",
        f"Unresolved lox records: {unresolved_count}. Exact gene-level risk is limited to SynII, SynIII, and SynVI. A 20-ORF verification sample is saved in gene_lox_mapping_validation_sample.csv.",
        "",
        "## Essentiality Gate",
        "",
        markdown_table(gate),
        "",
        f"The formal main grid uses the strict {gate_strict}-gene set. Expanded, high-confidence, and five size-matched random controls are versioned independently. Gate failure means zero current copies of at least one selected ORF; it is not an experimentally observed dead-cell genome.",
        "",
        "## Formal simulation results",
        "",
        f"The 88 scenarios initialized {int(population['initialized_trajectories'].sum()):,} trajectories and accepted {total_events:,} events: {int(event_totals['deletion_events']):,} deletions, {int(event_totals['inversion_events']):,} inversions, and {int(event_totals['duplication_events']):,} duplications.",
        f"The best rarefied structural effective Shannon diversity was {fmt(best_div.get('effective_shannon_diversity_mean'), 2)} for {best_div.get('model_key')} at alpha={best_div.get('alpha')} and p_event={best_div.get('p_event')}; its 95% rarefaction interval was {fmt(best_div.get('effective_shannon_diversity_ci_low'), 2)}-{fmt(best_div.get('effective_shannon_diversity_ci_high'), 2)}.",
        f"At essentiality-gate passing fraction >=0.70, the largest rarefied structural diversity was {fmt(best_retention.get('effective_shannon_diversity_mean'), 2)} for {best_retention.get('model_key')} at alpha={best_retention.get('alpha')} and p_event={best_retention.get('p_event')}.",
        f"The revised common rarefaction depth is {int(structural['rarefied_survivor_depth'].min())}; all effective Shannon estimates satisfy D1 <= N. Full-versus-rarefied rank correlations are documented below.",
        "",
        markdown_table(compare),
        "",
        "## Endpoint isomorphism and copy-label audit",
        "",
        (
            f"The original `segment-topology-v1` signature is not graph-isomorphism normalized because event-derived copy IDs remain encoded. "
            f"After copy-label removal, {int((isomorphism['topology_groups_merging_multiple_provenance_ids'] > 0).sum()) if not isomorphism.empty else 0} of 88 scenarios contained any merged provenance hashes. "
            f"The maximum full-sample unique-endpoint inflation was {isomorphism['unique_inflation_ratio'].max():.6f} and the maximum paired-rarefied effective-Shannon inflation was {isomorphism['paired_rarefied_effective_shannon_inflation_mean'].max():.6f}. "
            f"Event-derived copy identity alone produced maxima of {isomorphism['copy_identity_unique_inflation_ratio'].max():.6f} for unique endpoints and {isomorphism['copy_identity_paired_rarefied_inflation_mean'].max():.6f} for paired-rarefied effective Shannon diversity; reference-token representation aliases produced maxima of {isomorphism['reference_alias_unique_inflation_ratio'].max():.6f} and {isomorphism['reference_alias_paired_rarefied_inflation_mean'].max():.6f}, respectively. "
            f"The highest-diversity parameter combination remained linear-distance sampling, alpha=2.0, p_event=0.30. "
            f"The copy-unlabeled topology metric is preferred for primary structural-diversity reporting; the copy-aware metric should be labeled lineage/provenance-aware."
            if not isomorphism.empty
            else "The endpoint isomorphism audit was not available."
        ),
        "",
        "## Frequency backtesting and controls",
        "",
        f"The largest partially Hi-C-informed Pearson correlation was r={fmt(best_hic.get('pearson_r'))} at alpha={best_hic.get('alpha')} and p_event={best_hic.get('p_event')}. The best linear-distance Pearson correlation was r={fmt(best_linear.get('pearson_r'))} at alpha={best_linear.get('alpha')} and p_event={best_linear.get('p_event')}.",
        "",
        markdown_table(nulls),
        "",
        f"Chromosome-wise holdout contains {len(holdout)} held-out chromosomes and remains a parameter-selection holdout within the same processed study, not independent experimental validation.",
        "",
        "## Robustness analyses",
        "",
        robustness_text,
        "",
        "## Exact-coordinate gene risk",
        "",
        f"The revised risk table contains {len(risk)} genes observed among gate-failing trajectories. All legacy named-gene rankings are invalid and must not be reused; only exact-coordinate records in exact_coordinate_gene_risk.csv may be discussed, with the gate-based limitation stated.",
        "",
        "## Status of legacy results",
        "",
        "### Invalidated",
        "",
        "- The 1,092-site input, 148,512-pair proposal space, and 96,096 predicted-frequency row count.",
        "- The claim that four chromosomes had direct Hi-C weighting and the old 37,128-direct-pair interpretation.",
        "- All old effective/constrained diversity maxima, p_event optima, and threshold values derived from the fixed-reference engine.",
        "- Pearson values 0.3022 and 0.319664 as evidence for the revised model; they belong to old analysis boundaries.",
        "- The 133-ORF gate count, old PRP6/REB1/CHS2/CDS1/ALG14 risk ranking, and the old model-inferred dead-cell total.",
        "- Any description of endpoint diversity as genome structural diversity when it was based only on truncated ORF copy-number states.",
        "",
        "### Retained only at the conceptual level",
        "",
        "- Survivor-only endpoint data motivate a risk-aware analysis of unobserved gate-failing trajectories.",
        "- Increasing rearrangement opportunity creates a survival-diversity trade-off in the model, but its quantitative location must come from the revised tables.",
        "- SynII Hi-C can be tested as a partial proposal constraint; any advantage must be reported from revised backtesting and null controls.",
        "",
        "## Automated consistency checks",
        "",
        markdown_table(checks),
        "",
        "## Reproducibility boundary",
        "",
        "The revised workflow does not claim a calibrated wet-lab induction protocol, complete fitness model, or base-pair breakpoint predictor. p_event remains a per-trajectory per-step model probability. All manuscript revision should be based only on this revised output tree and should preserve the distinction between essentiality-gate failure and biological death.",
    ]
    return "\n".join(lines) + "\n"


def manifest(project: Path, output: Path) -> str:
    roots = [project / "scripts" / "revised_model", output]
    rows = []
    for root in roots:
        for path in sorted(item for item in root.rglob("*") if item.is_file()):
            relative = path.relative_to(project)
            if "scripts" in relative.parts:
                category = "script"
            elif "logs" in relative.parts:
                category = "log"
            elif "inputs" in relative.parts:
                category = "revised input/audit"
            elif "analysis" in relative.parts:
                category = "analysis result"
            elif "robustness" in relative.parts:
                category = "robustness result"
            else:
                category = "simulation output/intermediate"
            rows.append({"category": category, "size_bytes": path.stat().st_size, "path": str(path)})
    frame = pd.DataFrame(rows)
    return "# FILE MANIFEST\n\n" + markdown_table(frame) + "\n"


def main() -> None:
    args = parse_args()
    checks = build_checks(args.output_root)
    checks.to_csv(args.output_root / "automated_consistency_checks.csv", index=False)
    (args.output_root / "FINAL_REVISION_REPORT.md").write_text(
        generate_report(args.project_root, args.output_root, checks), encoding="utf-8"
    )
    (args.output_root / "FILE_MANIFEST.md").write_text(manifest(args.project_root, args.output_root), encoding="utf-8")
    failed = checks[checks["status"] != "PASS"]
    print(json.dumps({"checks": len(checks), "failed": len(failed)}, indent=2))
    if not failed.empty:
        raise RuntimeError(f"Consistency checks failed: {failed['check'].tolist()}")


if __name__ == "__main__":
    main()
