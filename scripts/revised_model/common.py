from __future__ import annotations

import hashlib
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


LOX_LEFT_ARM = "ATAACTTCGTATA"
LOX_RIGHT_ARM = "TATACGAAGTTAT"
LOX_LENGTH = 34


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    name: str | None = None
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                records[name] = []
            elif name is not None:
                records[name].append(line.upper())
    return {key: "".join(parts) for key, parts in records.items()}


def locate_loxpsym(sequence: str) -> list[dict[str, Any]]:
    hits: list[dict[str, Any]] = []
    start = 0
    while True:
        index = sequence.find(LOX_LEFT_ARM, start)
        if index < 0:
            break
        candidate = sequence[index : index + LOX_LENGTH]
        if len(candidate) == LOX_LENGTH and candidate[-len(LOX_RIGHT_ARM) :] == LOX_RIGHT_ARM:
            hits.append(
                {
                    "position_bp": index + 1,
                    "end_bp": index + LOX_LENGTH,
                    "sequence": candidate,
                    "spacer": candidate[len(LOX_LEFT_ARM) : -len(LOX_RIGHT_ARM)],
                }
            )
        start = index + 1
    return hits


@dataclass(frozen=True)
class GenBankFeature:
    feature_type: str
    start: int
    end: int
    strand: str
    qualifiers: dict[str, list[str]]


@dataclass(frozen=True)
class GenBankRecord:
    accession: str
    sequence: str
    features: tuple[GenBankFeature, ...]


def _parse_qualifiers(lines: list[str]) -> dict[str, list[str]]:
    qualifiers: dict[str, list[str]] = {}
    key: str | None = None
    value_parts: list[str] = []

    def flush() -> None:
        nonlocal key, value_parts
        if key is None:
            return
        value = " ".join(value_parts).strip().strip('"')
        qualifiers.setdefault(key, []).append(value)
        key = None
        value_parts = []

    for raw in lines:
        text = raw.strip()
        if text.startswith("/") and "=" in text:
            flush()
            key, value = text[1:].split("=", 1)
            value_parts = [value]
            if value.endswith('"') and value.count('"') >= 2:
                flush()
        elif text.startswith("/"):
            flush()
            qualifiers.setdefault(text[1:], []).append("true")
        elif key is not None:
            value_parts.append(text)
            if text.endswith('"'):
                flush()
    flush()
    return qualifiers


def read_genbank(path: Path) -> GenBankRecord:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    accession = path.stem
    in_features = False
    in_origin = False
    sequence_parts: list[str] = []
    raw_features: list[tuple[str, str, list[str]]] = []
    current_type: str | None = None
    current_location = ""
    current_qualifiers: list[str] = []

    def flush_feature() -> None:
        nonlocal current_type, current_location, current_qualifiers
        if current_type is not None:
            raw_features.append((current_type, current_location, current_qualifiers.copy()))
        current_type = None
        current_location = ""
        current_qualifiers = []

    for line in lines:
        if line.startswith("ACCESSION"):
            tokens = line.split()
            if len(tokens) > 1:
                accession = tokens[1]
        if line.startswith("FEATURES"):
            in_features = True
            continue
        if line.startswith("ORIGIN"):
            flush_feature()
            in_features = False
            in_origin = True
            continue
        if line.startswith("//"):
            in_origin = False
            continue
        if in_origin:
            sequence_parts.append("".join(re.findall(r"[A-Za-z]+", line)).upper())
            continue
        if not in_features:
            continue
        feature_key = line[5:21].strip() if len(line) >= 21 else ""
        payload = line[21:].strip() if len(line) >= 21 else ""
        if feature_key:
            flush_feature()
            current_type = feature_key
            current_location = payload
        elif payload.startswith("/"):
            current_qualifiers.append(payload)
        elif current_qualifiers:
            current_qualifiers.append(payload)
        elif current_type is not None:
            current_location += payload

    sequence = "".join(sequence_parts)
    features: list[GenBankFeature] = []
    for feature_type, location, qualifier_lines in raw_features:
        numbers = [int(token) for token in re.findall(r"\d+", location)]
        if not numbers:
            continue
        features.append(
            GenBankFeature(
                feature_type=feature_type,
                start=min(numbers),
                end=max(numbers),
                strand="-" if "complement" in location else "+",
                qualifiers=_parse_qualifiers(qualifier_lines),
            )
        )
    return GenBankRecord(accession=accession, sequence=sequence, features=tuple(features))


def qualifier_first(feature: GenBankFeature, *keys: str) -> str:
    for key in keys:
        values = feature.qualifiers.get(key, [])
        if values and values[0]:
            return values[0]
    return ""


def qualifier_sgd_id(feature: GenBankFeature) -> str:
    for value in feature.qualifiers.get("db_xref", []):
        if value.startswith("SGD:"):
            return value.split(":", 1)[1]
    return ""


def unique_find(sequence: str, query: str) -> tuple[int | None, str]:
    first = sequence.find(query)
    if first < 0:
        return None, "not_found"
    if sequence.find(query, first + 1) >= 0:
        return None, "multiple_matches"
    return first, "unique_exact"


def map_feature_sequence(
    source_sequence: str,
    target_sequence: str,
    start: int,
    end: int,
    anchor_length: int = 80,
) -> tuple[int | None, int | None, str]:
    query = source_sequence[start - 1 : end]
    if not query:
        return None, None, "empty_source_sequence"
    index, status = unique_find(target_sequence, query)
    if index is not None:
        return index + 1, index + len(query), status
    if len(query) < anchor_length * 2:
        return None, None, status
    left = query[:anchor_length]
    right = query[-anchor_length:]
    left_index, left_status = unique_find(target_sequence, left)
    right_index, right_status = unique_find(target_sequence, right)
    if left_index is None or right_index is None or right_index < left_index:
        return None, None, f"anchor_failed:{left_status}:{right_status}"
    mapped_end = right_index + anchor_length
    if mapped_end - left_index <= 0:
        return None, None, "anchor_order_invalid"
    return left_index + 1, mapped_end, "unique_flanking_anchors"


def choose_nearest_bin(position_bp: int, bin_size: int = 10_000) -> int:
    return int(round(position_bp / bin_size) * bin_size)


def combination_count(n: int) -> int:
    return math.comb(n, 2) if n >= 2 else 0


def semicolon(values: Iterable[str]) -> str:
    return ";".join(sorted({str(value) for value in values if str(value)}))
