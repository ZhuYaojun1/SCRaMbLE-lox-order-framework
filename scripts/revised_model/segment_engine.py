from __future__ import annotations

import copy
import hashlib
import json
import math
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd


@dataclass
class SegmentInstance:
    template_id: str
    copy_id: str
    orientation: int
    length_bp: int | None
    orf_ids: tuple[str, ...] = ()
    essential_orf_ids: tuple[str, ...] = ()
    feature_ids: tuple[str, ...] = ()

    def clone(self, copy_id: str) -> "SegmentInstance":
        return SegmentInstance(
            template_id=self.template_id,
            copy_id=copy_id,
            orientation=self.orientation,
            length_bp=self.length_bp,
            orf_ids=self.orf_ids,
            essential_orf_ids=self.essential_orf_ids,
            feature_ids=self.feature_ids,
        )


@dataclass
class LoxInstance:
    template_id: str
    copy_id: str
    orientation: int = 1
    reference_position_bp: int | None = None
    coordinate_status: str = "resolved"

    def clone(self, copy_id: str) -> "LoxInstance":
        return LoxInstance(
            template_id=self.template_id,
            copy_id=copy_id,
            orientation=self.orientation,
            reference_position_bp=self.reference_position_bp,
            coordinate_status=self.coordinate_status,
        )


Part = SegmentInstance | LoxInstance


@dataclass
class EventOutcome:
    event_id: int
    event_type: str
    chromosome: str
    left_lox_copy_id: str
    right_lox_copy_id: str
    left_lox_template_id: str
    right_lox_template_id: str
    selected_lox_separation: int
    affected_segment_copy_ids: tuple[str, ...]
    affected_orf_ids: tuple[str, ...]
    affected_essential_orf_ids: tuple[str, ...]
    lox_count_before: int
    lox_count_after: int
    segment_count_before: int
    segment_count_after: int


class ChromosomeTopology:
    def __init__(self, chromosome: str, parts: list[Part]) -> None:
        self.chromosome = chromosome
        self.parts = parts
        self.validate()

    def validate(self) -> None:
        if len(self.parts) < 3 or len(self.parts) % 2 == 0:
            raise ValueError(f"{self.chromosome}: topology must contain segment/lox alternation.")
        for index, part in enumerate(self.parts):
            expected_segment = index % 2 == 0
            if expected_segment != isinstance(part, SegmentInstance):
                raise ValueError(f"{self.chromosome}: invalid part at {index}: {type(part).__name__}")
        segment_ids = [part.copy_id for part in self.parts if isinstance(part, SegmentInstance)]
        lox_ids = [part.copy_id for part in self.parts if isinstance(part, LoxInstance)]
        if len(segment_ids) != len(set(segment_ids)):
            raise ValueError(f"{self.chromosome}: duplicate segment copy ID.")
        if len(lox_ids) != len(set(lox_ids)):
            raise ValueError(f"{self.chromosome}: duplicate lox copy ID.")

    def segments(self) -> list[SegmentInstance]:
        return [part for part in self.parts if isinstance(part, SegmentInstance)]

    def lox(self) -> list[LoxInstance]:
        return [part for part in self.parts if isinstance(part, LoxInstance)]

    def lox_index_by_copy_id(self) -> dict[str, int]:
        return {lox.copy_id: index for index, lox in enumerate(self.lox())}

    def segment_order(self) -> list[str]:
        return [segment.copy_id for segment in self.segments()]

    def segment_orientations(self) -> list[int]:
        return [segment.orientation for segment in self.segments()]

    def junctions(self) -> list[str]:
        junctions: list[str] = []
        for part_index in range(1, len(self.parts), 2):
            left = self.parts[part_index - 1]
            lox = self.parts[part_index]
            right = self.parts[part_index + 1]
            assert isinstance(left, SegmentInstance)
            assert isinstance(lox, LoxInstance)
            assert isinstance(right, SegmentInstance)
            junctions.append(f"{left.copy_id}>{lox.copy_id}>{right.copy_id}")
        return junctions

    def canonical_payload(self) -> str:
        part_tokens: list[str] = []
        for part in self.parts:
            if isinstance(part, SegmentInstance):
                part_tokens.append(f"S({part.template_id},{part.copy_id},{part.orientation:+d})")
            else:
                part_tokens.append(f"L({part.template_id},{part.copy_id},{part.orientation:+d})")
        junctions = ",".join(self.junctions())
        active_lox = ",".join(lox.copy_id for lox in self.lox())
        return f"{self.chromosome}|PARTS:{','.join(part_tokens)}|JUNCTIONS:{junctions}|ACTIVE_LOX:{active_lox}"

    def apply_event(self, event_type: str, left_index: int, right_index: int, event_id: int) -> EventOutcome:
        lox_list = self.lox()
        if not (0 <= left_index < right_index < len(lox_list)):
            raise ValueError("Selected lox indices must be distinct and ordered.")
        left_lox = lox_list[left_index]
        right_lox = lox_list[right_index]
        left_part_index = left_index * 2 + 1
        right_part_index = right_index * 2 + 1
        affected_parts = self.parts[left_part_index + 1 : right_part_index]
        affected_segments = [part for part in affected_parts if isinstance(part, SegmentInstance)]
        affected_orfs = sorted({gene for segment in affected_segments for gene in segment.orf_ids})
        affected_essential = sorted(
            {gene for segment in affected_segments for gene in segment.essential_orf_ids}
        )
        lox_before = len(lox_list)
        segment_before = len(self.segments())

        if event_type == "deletion":
            # Keep the left recombination site as the single junction site and remove
            # the intervening DNA together with the right boundary copy.
            del self.parts[left_part_index + 1 : right_part_index + 1]
        elif event_type == "inversion":
            inverted: list[Part] = []
            for part in reversed(self.parts[left_part_index + 1 : right_part_index]):
                part_copy = copy.deepcopy(part)
                part_copy.orientation *= -1
                inverted.append(part_copy)
            self.parts[left_part_index + 1 : right_part_index] = inverted
        elif event_type == "duplication":
            # Tandemly copy the current interval plus its right boundary. The left
            # boundary is shared, while every copied segment/lox receives a new ID.
            duplicated: list[Part] = []
            for offset, part in enumerate(self.parts[left_part_index + 1 : right_part_index + 1]):
                if isinstance(part, SegmentInstance):
                    duplicated.append(part.clone(f"{part.template_id}@e{event_id}s{offset}"))
                else:
                    duplicated.append(part.clone(f"{part.template_id}@e{event_id}l{offset}"))
            self.parts[right_part_index + 1 : right_part_index + 1] = duplicated
        else:
            raise ValueError(f"Unsupported event type: {event_type}")
        self.validate()
        return EventOutcome(
            event_id=event_id,
            event_type=event_type,
            chromosome=self.chromosome,
            left_lox_copy_id=left_lox.copy_id,
            right_lox_copy_id=right_lox.copy_id,
            left_lox_template_id=left_lox.template_id,
            right_lox_template_id=right_lox.template_id,
            selected_lox_separation=right_index - left_index,
            affected_segment_copy_ids=tuple(segment.copy_id for segment in affected_segments),
            affected_orf_ids=tuple(affected_orfs),
            affected_essential_orf_ids=tuple(affected_essential),
            lox_count_before=lox_before,
            lox_count_after=len(self.lox()),
            segment_count_before=segment_before,
            segment_count_after=len(self.segments()),
        )


class GenomeState:
    def __init__(
        self,
        chromosomes: dict[str, ChromosomeTopology],
        reference_payloads: dict[str, str] | None = None,
        reference_orf_counts: Counter[str] | None = None,
    ) -> None:
        self.chromosomes = chromosomes
        self.event_counter = 0
        self.modified_chromosomes: set[str] = set()
        self.reference_payloads = reference_payloads or {
            name: topology.canonical_payload() for name, topology in chromosomes.items()
        }
        if reference_orf_counts is None:
            counts: Counter[str] = Counter()
            for topology in chromosomes.values():
                for segment in topology.segments():
                    counts.update(segment.orf_ids)
            self.reference_orf_counts = counts
        else:
            self.reference_orf_counts = reference_orf_counts

    def validate(self) -> None:
        for chromosome in self.chromosomes.values():
            chromosome.validate()

    def orf_copy_numbers(self) -> Counter[str]:
        counts: Counter[str] = Counter()
        for chromosome in self.chromosomes.values():
            for segment in chromosome.segments():
                counts.update(segment.orf_ids)
        return counts

    def essential_copy_numbers(self) -> Counter[str]:
        counts: Counter[str] = Counter()
        for chromosome in self.chromosomes.values():
            for segment in chromosome.segments():
                counts.update(segment.essential_orf_ids)
        return counts

    def gate_missing(self, gate_genes: set[str]) -> list[str]:
        counts = self.orf_copy_numbers()
        return sorted(gene for gene in gate_genes if counts.get(gene, 0) <= 0)

    def apply_event(self, event_type: str, chromosome: str, left_index: int, right_index: int) -> EventOutcome:
        self.event_counter += 1
        outcome = self.chromosomes[chromosome].apply_event(
            event_type, left_index, right_index, self.event_counter
        )
        self.modified_chromosomes.add(chromosome)
        return outcome

    def canonical_structure(self) -> str:
        tokens: list[str] = []
        for name in sorted(self.chromosomes):
            if name in self.modified_chromosomes:
                tokens.append(self.chromosomes[name].canonical_payload())
            else:
                reference_hash = hashlib.md5(self.reference_payloads[name].encode("utf-8")).hexdigest()
                tokens.append(f"{name}|REFERENCE:{reference_hash}")
        chromosomes = "||".join(tokens)
        return f"segment-topology-v1||{chromosomes}"

    def structural_signature(self) -> tuple[str, str]:
        canonical = self.canonical_structure()
        return hashlib.md5(canonical.encode("utf-8")).hexdigest(), canonical

    def orf_copy_number_signature(self) -> tuple[str, str]:
        current = self.orf_copy_numbers()
        keys = set(current) | set(self.reference_orf_counts)
        delta = sorted(
            (gene, current.get(gene, 0))
            for gene in keys
            if current.get(gene, 0) != self.reference_orf_counts.get(gene, 0)
        )
        reference_hash = hashlib.md5(
            json.dumps(sorted(self.reference_orf_counts.items()), separators=(",", ":")).encode("utf-8")
        ).hexdigest()
        canonical = json.dumps(
            {"schema": "orf-copy-number-v1", "reference": reference_hash, "delta": delta},
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
        )
        return hashlib.md5(canonical.encode("utf-8")).hexdigest(), canonical


def _split_ids(value: Any) -> tuple[str, ...]:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ()
    return tuple(token for token in str(value).split(";") if token and token != "nan")


def load_reference_genome(inputs_dir: Path) -> GenomeState:
    lox = pd.read_csv(inputs_dir / "lox_sites_validated.csv")
    segments = pd.read_csv(inputs_dir / "reference_segments.csv")
    chromosomes: dict[str, ChromosomeTopology] = {}
    for chromosome, segment_group in segments.groupby("chromosome", sort=True):
        segment_group = segment_group.sort_values("segment_order")
        lox_group = lox[lox["chromosome"] == chromosome].sort_values("lox_order")
        parts: list[Part] = []
        segment_records = list(segment_group.to_dict("records"))
        lox_records = list(lox_group.to_dict("records"))
        if len(segment_records) != len(lox_records) + 1:
            raise RuntimeError(f"{chromosome}: reference segment/lox count mismatch.")
        for index, segment in enumerate(segment_records):
            start = pd.to_numeric(pd.Series([segment.get("start_bp")]), errors="coerce").iloc[0]
            end = pd.to_numeric(pd.Series([segment.get("end_bp")]), errors="coerce").iloc[0]
            length = int(end - start + 1) if pd.notna(start) and pd.notna(end) and end >= start else None
            parts.append(
                SegmentInstance(
                    template_id=str(segment["segment_id"]),
                    copy_id=f"{segment['segment_id']}@0",
                    orientation=1,
                    length_bp=length,
                    orf_ids=_split_ids(segment.get("orf_ids")),
                    essential_orf_ids=_split_ids(segment.get("essential_orf_ids")),
                    feature_ids=_split_ids(segment.get("important_feature_ids")),
                )
            )
            if index < len(lox_records):
                lox_record = lox_records[index]
                position = pd.to_numeric(pd.Series([lox_record.get("position_bp")]), errors="coerce").iloc[0]
                parts.append(
                    LoxInstance(
                        template_id=str(lox_record["lox_id"]),
                        copy_id=f"{lox_record['lox_id']}@0",
                        orientation=1,
                        reference_position_bp=int(position) if pd.notna(position) else None,
                        coordinate_status=str(lox_record.get("coordinate_status", "unresolved")),
                    )
                )
        chromosomes[chromosome] = ChromosomeTopology(chromosome, parts)
    return GenomeState(chromosomes)


class PairSampler:
    def __init__(
        self,
        model_name: str,
        alpha: float,
        direct_weights: pd.DataFrame | None = None,
    ) -> None:
        self.model_name = model_name
        self.alpha = float(alpha)
        self.direct_pairs: list[tuple[str, str, str, float]] = []
        self.direct_lookup: dict[tuple[str, str, str], float] = {}
        if direct_weights is not None and not direct_weights.empty:
            for row in direct_weights.to_dict("records"):
                chromosome = str(row["chromosome"])
                first = f"{row['lox_id_1']}@0"
                second = f"{row['lox_id_2']}@0"
                value = float(row["direct_contact_weight"])
                key = (chromosome, *sorted((first, second)))
                self.direct_pairs.append((chromosome, first, second, value))
                self.direct_lookup[key] = value
        self.direct_total = sum(row[3] for row in self.direct_pairs)
        direct_array = np.fromiter((row[3] for row in self.direct_pairs), dtype=float)
        self.direct_probabilities = direct_array / direct_array.sum() if direct_array.size else direct_array
        self.direct_cdf = np.cumsum(self.direct_probabilities)
        if self.direct_cdf.size:
            self.direct_cdf[-1] = 1.0
        self._distance_cache: dict[tuple[int, float], tuple[np.ndarray, np.ndarray, float]] = {}

    def _distance_distribution(self, n: int) -> tuple[np.ndarray, np.ndarray, float]:
        key = (n, self.alpha)
        cached = self._distance_cache.get(key)
        if cached is not None:
            return cached
        distances = np.arange(1, n, dtype=int)
        weights = (n - distances) / ((1.0 + distances) ** self.alpha)
        total = float(weights.sum())
        probabilities = weights / total if total > 0 else np.array([], dtype=float)
        cached = distances, probabilities, total
        self._distance_cache[key] = cached
        return cached

    def _fallback_totals(self, state: GenomeState) -> tuple[list[str], np.ndarray, float]:
        chromosomes: list[str] = []
        totals: list[float] = []
        for name, topology in state.chromosomes.items():
            n = len(topology.lox())
            if n < 2:
                continue
            if self.model_name == "uniform_random":
                total = n * (n - 1) / 2
            else:
                total = self._distance_distribution(n)[2]
            chromosomes.append(name)
            totals.append(total)
        array = np.asarray(totals, dtype=float)
        return chromosomes, array, float(array.sum())

    def _sample_fallback(self, state: GenomeState, rng: np.random.Generator) -> tuple[str, int, int, float]:
        chromosomes, totals, total = self._fallback_totals(state)
        if total <= 0:
            raise RuntimeError("No eligible lox pair remains.")
        chromosome = str(rng.choice(chromosomes, p=totals / total))
        n = len(state.chromosomes[chromosome].lox())
        if self.model_name == "uniform_random":
            first, second = sorted(rng.choice(n, size=2, replace=False).tolist())
            return chromosome, int(first), int(second), 1.0
        distances, probabilities, _ = self._distance_distribution(n)
        distance = int(rng.choice(distances, p=probabilities))
        first = int(rng.integers(0, n - distance))
        second = first + distance
        weight = 1.0 / ((1.0 + distance) ** self.alpha)
        return chromosome, first, second, weight

    def sample(self, state: GenomeState, rng: np.random.Generator) -> tuple[str, int, int, str]:
        if self.model_name in {"uniform_random", "linear_distance"} or self.direct_total <= 0:
            chromosome, first, second, _ = self._sample_fallback(state, rng)
            return chromosome, first, second, "uniform" if self.model_name == "uniform_random" else "distance_fallback"

        _, _, fallback_total = self._fallback_totals(state)
        proposal_total = fallback_total + self.direct_total
        for _ in range(10_000):
            use_direct_component = rng.random() >= fallback_total / proposal_total
            if use_direct_component:
                selected = int(np.searchsorted(self.direct_cdf, rng.random(), side="right"))
                if selected >= len(self.direct_pairs):
                    selected = len(self.direct_pairs) - 1
                chromosome, first_copy, second_copy, contact = self.direct_pairs[selected]
                indices = state.chromosomes[chromosome].lox_index_by_copy_id()
                if first_copy not in indices or second_copy not in indices:
                    continue
                first, second = sorted((indices[first_copy], indices[second_copy]))
                return chromosome, first, second, "direct_hic"
            else:
                chromosome, first, second, _ = self._sample_fallback(state, rng)
                return chromosome, first, second, "distance_fallback"
        raise RuntimeError("Pair sampler exceeded rejection limit.")


def clone_reference(reference: GenomeState) -> GenomeState:
    # Reference parts are immutable during deletion/duplication; inversion deep-copies
    # every affected part before changing orientation. A shallow part-list clone avoids
    # copying roughly one thousand reference objects for every trajectory.
    chromosomes = {
        name: ChromosomeTopology(name, list(topology.parts))
        for name, topology in reference.chromosomes.items()
    }
    return GenomeState(
        chromosomes,
        reference_payloads=reference.reference_payloads,
        reference_orf_counts=reference.reference_orf_counts,
    )


def count_parts(state: GenomeState) -> tuple[int, int]:
    segments = sum(len(chromosome.segments()) for chromosome in state.chromosomes.values())
    lox = sum(len(chromosome.lox()) for chromosome in state.chromosomes.values())
    return segments, lox


def all_template_copy_counts(state: GenomeState) -> Counter[str]:
    counts: Counter[str] = Counter()
    for chromosome in state.chromosomes.values():
        counts.update(segment.template_id for segment in chromosome.segments())
    return counts
