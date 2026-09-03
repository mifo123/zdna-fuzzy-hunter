#!/usr/bin/env python3
"""Build Supplementary Table S1 from the released Shin analysis features."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

import zdna_fuzzy_detector as zfd  # noqa: E402


DEFAULT_INPUT = REPO_ROOT / "validation" / "data" / "shin_analysis_features.csv"
DEFAULT_OUTPUT = REPO_ROOT / "supplementary" / "Supplementary_Table_S1_Shin_per_locus.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def ranks(rows: list[dict[str, object]], field: str) -> dict[int, int]:
    order = sorted(
        range(len(rows)),
        key=lambda index: (-float(rows[index][field]), str(rows[index]["candidate_id"])),
    )
    return {index: rank for rank, index in enumerate(order, start=1)}


def main() -> None:
    args = parse_args()
    source_rows = read_csv(args.input.resolve())
    params = zfd.TunedParams()
    scored: list[dict[str, object]] = []
    for source in source_rows:
        features = zfd.features_from_table_row(source)
        score, components, active_rules = zfd.score_features(features, params)
        if float(source["grid_hit_count"]) <= 0.0:
            score = -1.0
        tss_distance = int(float(source["tss_distance_bp"]))
        overlap = float(source["grid_max_overlap_pct"])
        scored.append(
            {
                **source,
                "fuzzy_score": score,
                "components": components,
                "active_rules": active_rules,
                "balanced_prediction": int(
                    zfd.classify("balanced", score, components, tss_distance, overlap)
                    and score >= zfd.FUZZY_THRESHOLD
                ),
                "moderate_prediction": int(
                    zfd.classify("moderate", score, components, tss_distance, overlap)
                    and score >= zfd.FUZZY_THRESHOLD
                ),
                "strict_prediction": int(
                    zfd.classify("strict", score, components, tss_distance, overlap)
                    and score >= zfd.FUZZY_THRESHOLD
                ),
            }
        )

    fuzzy_rank = ranks(scored, "fuzzy_score")
    weighted_rank = ranks(scored, "grid_weighted_hit_fraction")
    overlap_rank = ranks(scored, "grid_max_overlap_pct")
    fields = [
        "candidate_id",
        "species",
        "chromosome",
        "start",
        "end",
        "shin_label",
        "fuzzy_score",
        "balanced_prediction",
        "moderate_prediction",
        "strict_prediction",
        "rank_fuzzy",
        "hunter_weighted_hit_fraction",
        "rank_hunter_weighted",
        "rank_shift_vs_hunter_weighted",
        "hunter_max_overlap_pct",
        "rank_hunter_overlap",
        "rank_shift_vs_hunter_overlap",
        "hunter_max_raw_score",
        "tss_distance_bp",
        "grid_hit_m2_l6_t30",
        "grid_overlap_pct_m2_l6_t30",
        "grid_hit_m2_l10_t60",
        "grid_overlap_pct_m2_l10_t60",
        "component_signal",
        "component_context",
        "component_evidence",
        "component_bias",
        "component_feasibility",
        "top_rules",
    ]
    output_rows: list[dict[str, object]] = []
    for index, row in enumerate(scored):
        components = row["components"]
        active_rules = row["active_rules"]
        output_rows.append(
            {
                "candidate_id": row["candidate_id"],
                "species": row["species"],
                "chromosome": row["chromosome"],
                "start": row["start"],
                "end": row["end"],
                "shin_label": row["shin_label"],
                "fuzzy_score": f"{float(row['fuzzy_score']):.4f}",
                "balanced_prediction": row["balanced_prediction"],
                "moderate_prediction": row["moderate_prediction"],
                "strict_prediction": row["strict_prediction"],
                "rank_fuzzy": fuzzy_rank[index],
                "hunter_weighted_hit_fraction": row["grid_weighted_hit_fraction"],
                "rank_hunter_weighted": weighted_rank[index],
                "rank_shift_vs_hunter_weighted": weighted_rank[index] - fuzzy_rank[index],
                "hunter_max_overlap_pct": row["grid_max_overlap_pct"],
                "rank_hunter_overlap": overlap_rank[index],
                "rank_shift_vs_hunter_overlap": overlap_rank[index] - fuzzy_rank[index],
                "hunter_max_raw_score": row["grid_max_raw_score"],
                "tss_distance_bp": row["tss_distance_bp"],
                "grid_hit_m2_l6_t30": row["grid_hit__m2_l6_t30"],
                "grid_overlap_pct_m2_l6_t30": row["grid_overlap_pct__m2_l6_t30"],
                "grid_hit_m2_l10_t60": row["grid_hit__m2_l10_t60"],
                "grid_overlap_pct_m2_l10_t60": row["grid_overlap_pct__m2_l10_t60"],
                "component_signal": f"{float(components['signal']):.4f}",
                "component_context": f"{float(components['context']):.4f}",
                "component_evidence": f"{float(components['evidence']):.4f}",
                "component_bias": f"{float(components['bias']):.4f}",
                "component_feasibility": f"{float(components['feasibility']):.4f}",
                "top_rules": ";".join(
                    f"{item['rule']}:{float(item['activation']):.3f}"
                    for item in active_rules[:3]
                ),
            }
        )

    output_path = args.output.resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"Wrote {len(output_rows)} rows to {output_path}")


if __name__ == "__main__":
    main()
