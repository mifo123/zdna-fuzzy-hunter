#!/usr/bin/env python3
"""Rebuild the analysis-ready Shin feature table from released inputs.

The public HG Shin FASTA files define the labelled source cohort.  The released
context table preserves the author-derived, non-Hunter variables used when the
fuzzy model was fitted.  Z-DNA Hunter features are recomputed locally with the
public implementation in this repository.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Dict


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

import zdna_fuzzy_detector as zfd  # noqa: E402


DATA_DIR = REPO_ROOT / "validation" / "data"
DEFAULT_CONTEXT = DATA_DIR / "shin_context_features.csv"
DEFAULT_SOURCE_DIR = DATA_DIR / "source"
DEFAULT_OUTPUT = DATA_DIR / "shin_analysis_features.csv"
SOURCE_FILES = {"positive": "ZFS_yes.fasta", "negative": "ZFS_no.fasta"}
CANONICAL_CHROMS = {f"chr{number}" for number in range(1, 23)} | {"chrX", "chrY"}
HEADER_PATTERN = re.compile(r"^(ZFS_\d+)_(.+)_(\d+)_(\d+)$")
CONTEXT_COLUMNS = [
    "candidate_id",
    "species",
    "length_bp",
    "evidence_type",
    "zdna_score",
    "tss_distance_bp",
    "promoter_overlap",
    "regulatory_marks",
    "repeat_overlap_pct",
    "primer_uniqueness",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--context", type=Path, default=DEFAULT_CONTEXT)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE_DIR)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--extract-context-from",
        type=Path,
        help="One-time maintenance action: extract non-grid columns from an existing feature table.",
    )
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(rows: list[dict[str, object]], path: Path, columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows([{column: row.get(column, "") for column in columns} for row in rows])


def extract_context(source: Path, output: Path) -> None:
    rows = read_csv(source)
    if not rows:
        raise SystemExit(f"No rows found in {source}")
    missing = [column for column in CONTEXT_COLUMNS if column not in rows[0]]
    if missing:
        raise SystemExit(f"Cannot extract required context columns: {missing}")
    write_csv(rows, output, CONTEXT_COLUMNS)
    print(f"Extracted {len(rows)} context rows to {output}")


def parse_source_header(header: str) -> tuple[str, str, int, int]:
    match = HEADER_PATTERN.fullmatch(header)
    if match is None:
        raise ValueError(f"Unexpected Shin FASTA header: {header}")
    fasta_id, chrom, start_text, end_text = match.groups()
    return fasta_id, chrom, int(start_text), int(end_text)


def source_records(source_dir: Path) -> Dict[str, dict[str, object]]:
    records: Dict[str, dict[str, object]] = {}
    for label, filename in SOURCE_FILES.items():
        for header, sequence in zfd.read_fasta(source_dir / filename):
            fasta_id, chrom, start_1based, end_1based = parse_source_header(header)
            if len(sequence) != end_1based - start_1based + 1:
                raise ValueError(
                    f"Sequence length does not match 1-based closed header coordinates: {header}"
                )
            candidate_id = f"SHIN_REAL_{label.upper()}_{fasta_id}"
            if candidate_id in records:
                raise ValueError(f"Duplicate source identifier: {candidate_id}")
            records[candidate_id] = {
                "candidate_id": candidate_id,
                "fasta_id": fasta_id,
                "label": label,
                "header": header,
                "chromosome": chrom,
                "start": start_1based - 1,
                "end": end_1based,
                "sequence": sequence.upper(),
            }
    return records


def format_float(value: float, digits: int) -> str:
    return f"{value:.{digits}f}"


def hunter_features(sequence: str) -> dict[str, object]:
    configs = zfd.HUNTER_PRESETS["shin-publication"]
    total_weight = sum(zfd.config_strictness(config) for config in configs)
    hit_count = 0
    strict_hit_count = 0
    weighted_hit = 0.0
    overlaps: list[float] = []
    scores: list[float] = []
    best_interval_length = 0
    strictest_threshold = 0.0
    longest_min_size = 0
    result: dict[str, object] = {}

    for config in configs:
        config_id = str(config["config_id"])
        hits = list(zfd.local_hunter_hits(sequence, config))
        hit = bool(hits)
        best_hit = max(
            hits,
            key=lambda item: (item.end - item.start, item.score_percent),
            default=None,
        )
        best_length = 0 if best_hit is None else best_hit.end - best_hit.start
        max_score = 0.0 if best_hit is None else best_hit.score_percent
        overlap_pct = 100.0 * best_length / max(1, len(sequence))
        result[f"grid_hit__{config_id}"] = int(hit)
        result[f"grid_overlap_pct__{config_id}"] = format_float(overlap_pct, 4)
        if not hit:
            continue
        hit_count += 1
        weighted_hit += zfd.config_strictness(config)
        overlaps.append(overlap_pct)
        scores.append(max_score)
        best_interval_length = max(best_interval_length, best_length)
        strictest_threshold = max(strictest_threshold, float(config["threshold"]))
        longest_min_size = max(longest_min_size, int(config["min_sequence_size"]))
        if float(config["threshold"]) >= 60.0 or int(config["min_sequence_size"]) >= 10:
            strict_hit_count += 1

    config_count = len(configs)
    hit_fraction = hit_count / config_count
    weighted_fraction = weighted_hit / total_weight
    strict_fraction = strict_hit_count / config_count
    if hit_count == 0:
        detection_class = "no_grid_hit"
    elif weighted_fraction >= 0.25:
        detection_class = "stable_parameter_hit"
    else:
        detection_class = "permissive_only"

    result.update(
        {
            "grid_config_count": config_count,
            "grid_hit_count": hit_count,
            "grid_hit_fraction": format_float(hit_fraction, 6),
            "grid_weighted_hit_fraction": format_float(weighted_fraction, 6),
            "grid_model1_hit_fraction": format_float(0.0, 6),
            "grid_model2_hit_fraction": format_float(hit_fraction, 6),
            "grid_model_consensus": format_float(0.0, 6),
            "grid_strict_hit_fraction": format_float(strict_fraction, 6),
            "grid_strictest_threshold": format_float(strictest_threshold, 2),
            "grid_longest_min_size": longest_min_size,
            "grid_max_overlap_pct": format_float(max(overlaps, default=0.0), 4),
            "grid_mean_overlap_pct": format_float(
                sum(overlaps) / len(overlaps) if overlaps else 0.0, 4
            ),
            "grid_max_raw_score": format_float(max(scores, default=0.0), 4),
            "grid_best_interval_length": best_interval_length,
            "grid_detection_class": detection_class,
        }
    )
    return result


def build(context_path: Path, source_dir: Path, output_path: Path) -> None:
    context_rows = read_csv(context_path)
    source = source_records(source_dir)
    included_source = {
        candidate_id: row
        for candidate_id, row in source.items()
        if row["chromosome"] in CANONICAL_CHROMS
    }
    excluded_source = [row for row in source.values() if row["chromosome"] not in CANONICAL_CHROMS]
    context_ids = {row["candidate_id"] for row in context_rows}
    if context_ids != set(included_source):
        missing = sorted(set(included_source) - context_ids)
        extra = sorted(context_ids - set(included_source))
        raise RuntimeError(f"Context/source cohort mismatch; missing={missing}, extra={extra}")

    output_rows: list[dict[str, object]] = []
    for context in context_rows:
        raw = source[context["candidate_id"]]
        start = int(raw["start"])
        end = int(raw["end"])
        length = int(float(context["length_bp"]))
        row: dict[str, object] = dict(context)
        row.update(
            {
                "chromosome": raw["chromosome"],
                "start": start,
                "end": start + length,
                "shin_label": raw["label"],
                "context_chromosome": raw["chromosome"],
                "context_window_start": start,
                "context_window_end": end,
                "shin_window_start": start,
                "shin_window_end": end,
                "shin_window_length": len(str(raw["sequence"])),
            }
        )
        row.update(hunter_features(str(raw["sequence"])))
        output_rows.append(row)

    analysis_columns = [
        "candidate_id",
        "species",
        "chromosome",
        "start",
        "end",
        "length_bp",
        "evidence_type",
        "zdna_score",
        "tss_distance_bp",
        "promoter_overlap",
        "regulatory_marks",
        "repeat_overlap_pct",
        "primer_uniqueness",
        "shin_label",
        "context_chromosome",
        "context_window_start",
        "context_window_end",
        "shin_window_start",
        "shin_window_end",
        "shin_window_length",
    ]
    grid_columns = [
        "grid_hit__m2_l6_t30",
        "grid_overlap_pct__m2_l6_t30",
        "grid_hit__m2_l10_t60",
        "grid_overlap_pct__m2_l10_t60",
        "grid_config_count",
        "grid_hit_count",
        "grid_hit_fraction",
        "grid_weighted_hit_fraction",
        "grid_model1_hit_fraction",
        "grid_model2_hit_fraction",
        "grid_model_consensus",
        "grid_strict_hit_fraction",
        "grid_strictest_threshold",
        "grid_longest_min_size",
        "grid_max_overlap_pct",
        "grid_mean_overlap_pct",
        "grid_max_raw_score",
        "grid_best_interval_length",
        "grid_detection_class",
    ]
    write_csv(output_rows, output_path, analysis_columns + grid_columns)
    positives = sum(row["label"] == "positive" for row in source.values())
    negatives = len(source) - positives
    print(
        f"Verified source cohort: {len(source)} loci "
        f"({positives} positive, {negatives} negative)."
    )
    print(
        f"Included {len(output_rows)} primary-chromosome loci; "
        f"excluded {len(excluded_source)} *_random-contig loci."
    )
    print(f"Wrote analysis feature table to {output_path}")


def main() -> None:
    args = parse_args()
    if args.extract_context_from is not None:
        extract_context(args.extract_context_from.resolve(), args.context.resolve())
        return
    build(args.context.resolve(), args.source_dir.resolve(), args.output.resolve())


if __name__ == "__main__":
    main()
