from __future__ import annotations

import argparse
import json
import math
import shutil
from bisect import bisect_left
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from common import (
    choose_nearest_bin,
    combination_count,
    locate_loxpsym,
    map_feature_sequence,
    qualifier_first,
    qualifier_sgd_id,
    read_fasta,
    read_genbank,
    semicolon,
    sha256_file,
)


CHROM_CONFIG = {
    "SynII": {"fasta": "chr2", "genbank": "CP013608.gb"},
    "SynIII": {"fasta": "chr3", "genbank": "KC880027.gb"},
    "SynVI": {"fasta": "chr6", "genbank": "CP135953.gb"},
    "SynIXR": {"fasta": "chr9", "genbank": "JN020955.gb"},
}
EXPECTED_VALID = {"SynII": 273, "SynIII": 100, "SynVI": 71, "SynIXR": 44}
REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build audited inputs for the revised SCRaMbLE model.")
    parser.add_argument("--project-root", type=Path, default=REPOSITORY_ROOT)
    parser.add_argument("--output-root", type=Path, default=REPOSITORY_ROOT / "outputs" / "revised_model")
    parser.add_argument("--genbank-source", type=Path, default=REPOSITORY_ROOT / "outputs" / "revised_model" / "input_sources")
    return parser.parse_args()


def parse_fig2c(workbook: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, list[int]]]:
    raw = pd.read_excel(workbook, sheet_name="Fig2C", header=None)
    orders = pd.to_numeric(raw.iloc[0, 2:], errors="coerce")
    records: list[dict[str, Any]] = []
    audits: list[dict[str, Any]] = []
    removed_examples: dict[str, list[int]] = {}
    for row_index in range(1, len(raw)):
        chromosome = str(raw.iat[row_index, 1]).strip()
        if chromosome not in EXPECTED_VALID:
            continue
        values = pd.to_numeric(raw.iloc[row_index, 2:], errors="coerce")
        valid = orders.notna() & values.notna()
        removed = orders[orders.notna() & values.isna()].astype(int).tolist()
        removed_examples[chromosome] = removed[:10]
        total = float(values[valid].sum())
        for order_value, count_value in zip(orders[valid], values[valid]):
            order = int(order_value)
            count = float(count_value)
            records.append(
                {
                    "chromosome": chromosome,
                    "lox_id": f"{chromosome}_loxPsym_{order}",
                    "lox_order": order,
                    "observed_rearrangement_count": count,
                    "observed_rearrangement_frequency": count / total if total > 0 else np.nan,
                    "source_file": str(workbook),
                    "source_sheet": "Fig2C",
                }
            )
        n_valid = int(valid.sum())
        audits.append(
            {
                "chromosome": chromosome,
                "raw_header_columns": int(orders.notna().sum()),
                "nonempty_records": n_valid,
                "removed_blank_records": int((orders.notna() & values.isna()).sum()),
                "unique_lox_id": n_valid,
                "duplicate_lox_id": 0,
                "lox_order_min": int(orders[valid].min()),
                "lox_order_max": int(orders[valid].max()),
                "candidate_same_chromosome_pairs": combination_count(n_valid),
                "sites_with_observed_frequency": n_valid,
                "removed_blank_examples": ";".join(map(str, removed[:10])),
            }
        )
    frame = pd.DataFrame(records).sort_values(["chromosome", "lox_order"]).reset_index(drop=True)
    audit = pd.DataFrame(audits).sort_values("chromosome").reset_index(drop=True)
    actual = frame.groupby("chromosome").size().to_dict()
    if actual != EXPECTED_VALID:
        raise RuntimeError(f"Fig2C nonempty record assertion failed: {actual} != {EXPECTED_VALID}")
    if frame["lox_id"].duplicated().any() or len(frame) != 488:
        raise RuntimeError("Validated Fig2C input must contain 488 unique lox IDs.")
    if int(audit["removed_blank_records"].sum()) != 604:
        raise RuntimeError("Expected exactly 604 padded blank records to be removed.")
    if int(audit["candidate_same_chromosome_pairs"].sum()) != 45_509:
        raise RuntimeError("Expected 45,509 same-chromosome candidate pairs.")
    return frame, audit, removed_examples


def attach_fasta_coordinates(lox: pd.DataFrame, fasta_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    sequences = read_fasta(fasta_path)
    rows: list[pd.DataFrame] = []
    evidence: list[dict[str, Any]] = []
    for chromosome, config in CHROM_CONFIG.items():
        chrom_rows = lox[lox["chromosome"] == chromosome].copy().sort_values("lox_order")
        sequence = sequences[config["fasta"]]
        hits = locate_loxpsym(sequence)
        evidence.append(
            {
                "chromosome": chromosome,
                "fasta_record": config["fasta"],
                "fasta_length": len(sequence),
                "fig2c_nonempty_records": len(chrom_rows),
                "physical_loxpsym_hits": len(hits),
                "coordinate_complete": len(hits) == len(chrom_rows),
                "unresolved_records": max(0, len(chrom_rows) - len(hits)),
            }
        )
        positions = [hit["position_bp"] for hit in hits]
        ends = [hit["end_bp"] for hit in hits]
        spacers = [hit["spacer"] for hit in hits]
        chrom_rows["position_bp"] = np.nan
        chrom_rows["end_bp"] = np.nan
        chrom_rows["lox_sequence_spacer"] = ""
        n_assign = min(len(chrom_rows), len(hits))
        assign_indices = chrom_rows.index[:n_assign]
        chrom_rows.loc[assign_indices, "position_bp"] = positions[:n_assign]
        chrom_rows.loc[assign_indices, "end_bp"] = ends[:n_assign]
        chrom_rows.loc[assign_indices, "lox_sequence_spacer"] = spacers[:n_assign]
        chrom_rows["coordinate_status"] = np.where(
            chrom_rows["position_bp"].notna(), "resolved_from_local_fasta_motif", "unresolved_no_matching_fasta_loxpsym"
        )
        chrom_rows["active_for_topology"] = True
        chrom_rows["eligible_for_bp_mapping"] = chrom_rows["position_bp"].notna()
        chrom_rows["fasta_record"] = config["fasta"]
        rows.append(chrom_rows)
    return pd.concat(rows).sort_values(["chromosome", "lox_order"]).reset_index(drop=True), pd.DataFrame(evidence)


def build_hic_coverage(lox: pd.DataFrame, hic_path: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    hic = pd.read_csv(hic_path)
    hic = hic[
        hic["chromosome_1"].astype(str).str.lower().eq("synii")
        & hic["chromosome_2"].astype(str).str.lower().eq("synii")
    ].copy()
    lookup: dict[tuple[int, int], float] = {}
    for row in hic.itertuples(index=False):
        low, high = sorted((int(row.position_1), int(row.position_2)))
        lookup[(low, high)] = float(row.contact_probability)

    coverage: list[dict[str, Any]] = []
    direct_rows: list[dict[str, Any]] = []
    scale_rows: list[dict[str, Any]] = []
    for chromosome, group in lox.groupby("chromosome", sort=True):
        records = list(group.sort_values("lox_order").to_dict("records"))
        total = combination_count(len(records))
        direct = zero = missing = unresolved = 0
        fallback_weights: list[float] = []
        direct_weights: list[float] = []
        for i in range(len(records) - 1):
            for j in range(i + 1, len(records)):
                first, second = records[i], records[j]
                fallback = 1.0 / (1.0 + (j - i))
                fallback_weights.append(fallback)
                if not (math.isfinite(float(first.get("position_bp", np.nan))) and math.isfinite(float(second.get("position_bp", np.nan)))):
                    unresolved += 1
                    continue
                if chromosome != "SynII":
                    missing += 1
                    continue
                bin_1 = choose_nearest_bin(int(first["position_bp"]))
                bin_2 = choose_nearest_bin(int(second["position_bp"]))
                key = tuple(sorted((bin_1, bin_2)))
                value = lookup.get(key)
                if value is None:
                    missing += 1
                elif not math.isfinite(value) or value <= 0:
                    zero += 1
                else:
                    direct += 1
                    direct_weights.append(value)
                    direct_rows.append(
                        {
                            "chromosome": chromosome,
                            "lox_id_1": first["lox_id"],
                            "lox_id_2": second["lox_id"],
                            "hic_bin_1": bin_1,
                            "hic_bin_2": bin_2,
                            "direct_contact_weight": value,
                            "source_file": str(hic_path),
                        }
                    )
        fallback = total - direct
        coverage.append(
            {
                "chromosome": chromosome,
                "total_candidate_pairs": total,
                "direct_hic_pairs": direct,
                "distance_fallback_pairs": fallback,
                "missing_contact_pairs": missing,
                "zero_contact_pairs": zero,
                "unresolved_coordinate_pairs": unresolved,
                "direct_hic_fraction": direct / total if total else 0.0,
                "fallback_fraction": fallback / total if total else 0.0,
                "model_label": "partially Hi-C-informed sampling with distance-based fallback",
            }
        )
        scale_rows.append(
            {
                "chromosome": chromosome,
                "direct_weight_min": min(direct_weights) if direct_weights else np.nan,
                "direct_weight_median": float(np.median(direct_weights)) if direct_weights else np.nan,
                "direct_weight_max": max(direct_weights) if direct_weights else np.nan,
                "fallback_alpha1_min": min(fallback_weights) if fallback_weights else np.nan,
                "fallback_alpha1_median": float(np.median(fallback_weights)) if fallback_weights else np.nan,
                "fallback_alpha1_max": max(fallback_weights) if fallback_weights else np.nan,
            }
        )
    coverage_frame = pd.DataFrame(coverage)
    totals = {
        column: int(coverage_frame[column].sum())
        for column in [
            "total_candidate_pairs",
            "direct_hic_pairs",
            "distance_fallback_pairs",
            "missing_contact_pairs",
            "zero_contact_pairs",
            "unresolved_coordinate_pairs",
        ]
    }
    coverage_frame = pd.concat(
        [
            coverage_frame,
            pd.DataFrame(
                [
                    {
                        "chromosome": "ALL",
                        **totals,
                        "direct_hic_fraction": totals["direct_hic_pairs"] / totals["total_candidate_pairs"],
                        "fallback_fraction": totals["distance_fallback_pairs"] / totals["total_candidate_pairs"],
                        "model_label": "partially Hi-C-informed sampling with distance-based fallback",
                    }
                ]
            ),
        ],
        ignore_index=True,
    )
    return coverage_frame, pd.DataFrame(direct_rows), pd.DataFrame(scale_rows)


def _sgd_maps(annotation_path: Path) -> tuple[dict[str, dict[str, str]], dict[str, dict[str, str]]]:
    annotation = pd.read_csv(annotation_path, low_memory=False)
    by_id: dict[str, dict[str, str]] = {}
    by_name: dict[str, dict[str, str]] = {}
    for row in annotation.to_dict("records"):
        record = {key: str(value) for key, value in row.items()}
        gene_id = record.get("gene_id", "")
        gene_name = record.get("gene_name", "")
        if gene_id and gene_id != "nan":
            by_id[gene_id] = record
        if gene_name and gene_name != "nan":
            by_name[gene_name.upper()] = record
    return by_id, by_name


def map_features(
    lox: pd.DataFrame,
    fasta_path: Path,
    genbank_dir: Path,
    annotation_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sequences = read_fasta(fasta_path)
    by_id, by_name = _sgd_maps(annotation_path)
    all_features: list[dict[str, Any]] = []
    source_audit: list[dict[str, Any]] = []
    for chromosome, config in CHROM_CONFIG.items():
        gb_path = genbank_dir / config["genbank"]
        record = read_genbank(gb_path)
        target = sequences[config["fasta"]]
        selected = [
            feature
            for feature in record.features
            if feature.feature_type in {"gene", "centromere", "rep_origin", "telomere"}
        ]
        mapped_count = 0
        for number, feature in enumerate(selected, start=1):
            mapped_start, mapped_end, status = map_feature_sequence(
                record.sequence, target, feature.start, feature.end
            )
            gene_name = qualifier_first(feature, "gene", "note", "locus_tag")
            sgd_id = qualifier_sgd_id(feature)
            annotation = by_id.get(sgd_id) or by_name.get(gene_name.upper()) or {}
            if mapped_start is not None:
                mapped_count += 1
            all_features.append(
                {
                    "synthetic_chromosome": chromosome,
                    "feature_record_id": f"{chromosome}_{feature.feature_type}_{number:04d}",
                    "feature_type": "ORF" if feature.feature_type == "gene" else feature.feature_type,
                    "sgd_gene_id": sgd_id or annotation.get("gene_id", ""),
                    "gene_name": gene_name or annotation.get("gene_name", ""),
                    "source_start_bp": feature.start,
                    "source_end_bp": feature.end,
                    "synthetic_start_bp": mapped_start,
                    "synthetic_end_bp": mapped_end,
                    "synthetic_midpoint_bp": (mapped_start + mapped_end) / 2 if mapped_start is not None else np.nan,
                    "strand": feature.strand,
                    "essential_status": annotation.get("essential_status", "Unknown"),
                    "mapping_status": status,
                    "source_accession": record.accession,
                    "source_file": str(gb_path),
                }
            )
        source_audit.append(
            {
                "chromosome": chromosome,
                "source_accession": record.accession,
                "source_sequence_length": len(record.sequence),
                "target_fasta_record": config["fasta"],
                "target_sequence_length": len(target),
                "selected_features": len(selected),
                "mapped_features": mapped_count,
                "unresolved_features": len(selected) - mapped_count,
            }
        )

    features = pd.DataFrame(all_features)
    mapped_rows: list[dict[str, Any]] = []
    for row in features.to_dict("records"):
        chromosome = row["synthetic_chromosome"]
        chrom_lox = lox[
            (lox["chromosome"] == chromosome) & lox["position_bp"].notna()
        ].sort_values("position_bp")
        complete = int(chrom_lox.shape[0]) == EXPECTED_VALID[chromosome]
        midpoint = row.get("synthetic_midpoint_bp")
        left_id = right_id = segment_id = ""
        risk_eligible = False
        interval_status = "feature_coordinate_unresolved"
        if midpoint is not None and pd.notna(midpoint):
            positions = chrom_lox["position_bp"].astype(int).tolist()
            insert = bisect_left(positions, float(midpoint))
            left_id = chrom_lox.iloc[insert - 1]["lox_id"] if insert > 0 else "TEL_LEFT"
            right_id = chrom_lox.iloc[insert]["lox_id"] if insert < len(chrom_lox) else "TEL_RIGHT"
            segment_id = f"{chromosome}_SEG_{insert:04d}"
            interval_status = "exact_adjacent_lox" if complete else "chromosome_lox_coordinate_incomplete"
            risk_eligible = complete and row["feature_type"] == "ORF"
        mapped_rows.append(
            {
                **row,
                "left_lox_id": left_id,
                "right_lox_id": right_id,
                "segment_id": segment_id,
                "interval_mapping_status": interval_status,
                "eligible_for_gene_risk": risk_eligible,
            }
        )
    feature_mapping = pd.DataFrame(mapped_rows)
    genes = feature_mapping[feature_mapping["feature_type"] == "ORF"].copy()
    genes = genes[
        [
            "synthetic_chromosome",
            "sgd_gene_id",
            "gene_name",
            "synthetic_start_bp",
            "synthetic_end_bp",
            "synthetic_midpoint_bp",
            "strand",
            "essential_status",
            "left_lox_id",
            "right_lox_id",
            "segment_id",
            "mapping_status",
            "interval_mapping_status",
            "eligible_for_gene_risk",
            "source_accession",
            "source_file",
        ]
    ].sort_values(["synthetic_chromosome", "synthetic_start_bp"], na_position="last")
    eligible = genes[genes["eligible_for_gene_risk"]].copy()
    sample_n = min(20, len(eligible))
    sample = eligible.sample(sample_n, random_state=20260717).sort_values(
        ["synthetic_chromosome", "synthetic_start_bp"]
    )
    return genes.reset_index(drop=True), feature_mapping.reset_index(drop=True), pd.DataFrame(source_audit), sample.reset_index(drop=True)


def build_reference_segments(lox: pd.DataFrame, feature_mapping: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for chromosome, chrom_lox in lox.groupby("chromosome", sort=True):
        ordered = chrom_lox.sort_values("lox_order").to_dict("records")
        for segment_index in range(len(ordered) + 1):
            left = "TEL_LEFT" if segment_index == 0 else ordered[segment_index - 1]["lox_id"]
            right = "TEL_RIGHT" if segment_index == len(ordered) else ordered[segment_index]["lox_id"]
            left_bp = 1 if segment_index == 0 else ordered[segment_index - 1].get("end_bp")
            right_bp = np.nan if segment_index == len(ordered) else ordered[segment_index].get("position_bp")
            segment_id = f"{chromosome}_SEG_{segment_index:04d}"
            features = feature_mapping[feature_mapping["segment_id"] == segment_id]
            rows.append(
                {
                    "chromosome": chromosome,
                    "segment_id": segment_id,
                    "segment_order": segment_index,
                    "left_lox_id": left,
                    "right_lox_id": right,
                    "start_bp": left_bp,
                    "end_bp": right_bp,
                    "coordinate_status": "resolved" if pd.notna(left_bp) and pd.notna(right_bp) else "partly_unresolved",
                    "orf_ids": semicolon(features.loc[features["feature_type"] == "ORF", "sgd_gene_id"].astype(str)),
                    "orf_names": semicolon(features.loc[features["feature_type"] == "ORF", "gene_name"].astype(str)),
                    "essential_orf_ids": semicolon(
                        features.loc[
                            (features["feature_type"] == "ORF")
                            & (features["essential_status"] == "Essential")
                            & (features["eligible_for_gene_risk"]),
                            "sgd_gene_id",
                        ].astype(str)
                    ),
                    "important_feature_ids": semicolon(
                        features.loc[features["feature_type"] != "ORF", "feature_record_id"].astype(str)
                    ),
                }
            )
    return pd.DataFrame(rows)


def write_audit_markdown(audit: pd.DataFrame, evidence: pd.DataFrame, path: Path) -> None:
    def markdown_table(frame: pd.DataFrame) -> str:
        columns = [str(column) for column in frame.columns]
        lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join(["---"] * len(columns)) + " |"]
        for row in frame.fillna("").astype(str).itertuples(index=False, name=None):
            lines.append("| " + " | ".join(value.replace("|", "\\|") for value in row) + " |")
        return "\n".join(lines)

    lines = [
        "# Revised lox input audit",
        "",
        "The Zhang 2022 Fig. 2C parser now requires both a valid header order and a non-empty chromosome-specific value.",
        "The 604 padded blank cells are excluded from the active topology, candidate-pair construction, and predicted-frequency output.",
        "",
        "## Fig. 2C row audit",
        "",
        markdown_table(audit),
        "",
        "## Physical-coordinate evidence",
        "",
        markdown_table(evidence),
        "",
        "The supplied local FASTA resolves 487 of 488 Fig. 2C records. SynIXR contains 44 non-empty experimental records but 43 verifiable loxPsym motifs; the unmatched record is retained as an order-defined topology node with unresolved bp coordinate and is excluded from gene-level risk mapping.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_readme(output_root: Path, provenance: dict[str, Any]) -> None:
    text = f"""# Revised SCRaMbLE model inputs

This directory is independent of the original processed inputs and results. Original project files were read only.

## Lox input

Fig. 2C cells are retained only when both the shared lox-order header and the chromosome-specific event-count cell are non-empty. This yields 273 SynII, 100 SynIII, 71 SynVI, and 44 SynIXR records (488 total; 45,509 same-chromosome pairs). The local Zhou FASTA provides motif-supported bp coordinates for 487 records. One SynIXR record remains unresolved and is not interpolated.

## Hi-C matching

Only the synII matrix in Fig. 5A is used. ring_synII is a separate matrix and is not silently assigned to another chromosome. Exact SynII lox bp coordinates are rounded to the nearest 10-kb matrix bin. A positive finite matrix entry is a direct weight; all other pairs use the current-structure distance fallback. The model label is **partially Hi-C-informed sampling with distance-based fallback**.

Direct contacts and fallback weights both lie in (0, 1]. The fallback is 1/(1+d)^alpha, where d is current topological lox separation, so its scale remains defined after deletion, inversion, and duplication. For an eligible pair, the target proposal mass is the fallback weight plus the direct contact weight when that exact SynII reference pair is available. The sampler draws these additive components directly; deleted reference pairs are rejected from the direct component and the remaining eligible mass is renormalized. Duplicated lox copies and all non-SynII pairs receive only the distance component. Weight ranges are reported in `hic_weight_scale_audit.csv`.

## Gene and feature mapping

Synthetic-chromosome features are taken from accessioned GenBank records and mapped to the supplied local FASTA by unique full-sequence matching or unique flanking anchors. Gene intervals are assigned to adjacent motif-supported lox coordinates. No proportional chromosome-wide scaling is used. SynIXR gene-risk eligibility is disabled because its Fig. 2C record count and physical motif count disagree.

## Provenance

```json
{json.dumps(provenance, indent=2, sort_keys=True)}
```
"""
    (output_root / "README.md").write_text(text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    output = args.output_root
    inputs = output / "inputs"
    sources = output / "input_sources"
    inputs.mkdir(parents=True, exist_ok=True)
    sources.mkdir(parents=True, exist_ok=True)

    workbook = args.project_root / "data" / "raw" / "zhang_2022_gce" / "41467_2022_33606_MOESM4_ESM.xlsx"
    fasta = next((args.project_root / "data" / "raw" / "zhou_2023_nsr").rglob("GSE168182_yZSJ025.fa"))
    hic = args.project_root / "data" / "processed" / "real_hic_contact_matrix.csv"
    annotation = args.project_root / "data" / "processed" / "real_gene_annotation.csv"

    source_files = [workbook, fasta, hic, annotation]
    for name in ["CP013608.gb", "KC880027.gb", "CP135953.gb", "JN020955.gb"]:
        source = args.genbank_source / name
        if not source.exists():
            raise FileNotFoundError(f"Missing accessioned GenBank source: {source}")
        destination = sources / name
        shutil.copy2(source, destination)
        source_files.append(destination)

    lox, audit, _ = parse_fig2c(workbook)
    lox, coordinate_evidence = attach_fasta_coordinates(lox, fasta)
    lox.to_csv(inputs / "lox_sites_validated.csv", index=False)
    audit.to_csv(inputs / "lox_input_audit.csv", index=False)
    coordinate_evidence.to_csv(inputs / "lox_coordinate_evidence.csv", index=False)
    write_audit_markdown(audit, coordinate_evidence, inputs / "lox_input_audit.md")

    coverage, direct_weights, scale_audit = build_hic_coverage(lox, hic)
    coverage.to_csv(inputs / "hic_pair_coverage.csv", index=False)
    direct_weights.to_csv(inputs / "hic_direct_pair_weights.csv", index=False)
    scale_audit.to_csv(inputs / "hic_weight_scale_audit.csv", index=False)

    genes, features, source_audit, validation_sample = map_features(
        lox, fasta, sources, annotation
    )
    genes.to_csv(inputs / "gene_lox_exact_mapping.csv", index=False)
    features.to_csv(inputs / "feature_lox_exact_mapping.csv", index=False)
    source_audit.to_csv(inputs / "genbank_mapping_audit.csv", index=False)
    validation_sample.to_csv(inputs / "gene_lox_mapping_validation_sample.csv", index=False)
    segments = build_reference_segments(lox, features)
    segments.to_csv(inputs / "reference_segments.csv", index=False)

    provenance = {
        "builder": str(Path(__file__).resolve()),
        "sources": {str(path): sha256_file(path) for path in source_files},
        "assertions": {
            "validated_lox_records": len(lox),
            "removed_blank_records": int(audit["removed_blank_records"].sum()),
            "candidate_same_chromosome_pairs": int(audit["candidate_same_chromosome_pairs"].sum()),
            "coordinate_resolved_lox": int(lox["position_bp"].notna().sum()),
            "gene_risk_eligible_orfs": int(genes["eligible_for_gene_risk"].sum()),
        },
    }
    (inputs / "input_provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    write_readme(output, provenance)
    print(json.dumps(provenance["assertions"], indent=2))


if __name__ == "__main__":
    main()
