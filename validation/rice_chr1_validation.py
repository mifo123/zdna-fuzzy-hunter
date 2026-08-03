#!/usr/bin/env python3
"""Frozen-model check against the original rice GSE252954 ZIP-seq tracks.

The primary endpoint is defined before signal inspection. Predictions are
sampled deterministically across fuzzy-score deciles. A 500-bp window is
experimentally positive only when both ZIP-seq IP replicates exceed both the
input and IgG controls and mean IP is in the top 10% of sampled windows.
CUT&Tag is retained as a secondary continuous, orthogonal signal check. No
GSE252954 value is used to tune model weights, rules or the decision threshold.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import statistics
import sys
from pathlib import Path
from typing import Sequence

from u2os_chr22_validation import (
    bootstrap_rate_difference,
    mean_signal,
    pearson,
    percentile,
    rate,
    spearman,
    write_csv,
)


SEED = 20260805
SAMPLE_PER_DECILE = 2000
WINDOW_BP = 500
PSEUDOCOUNT = 1.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--igg", type=Path, required=True)
    parser.add_argument("--ip-rep1", type=Path, required=True)
    parser.add_argument("--ip-rep2", type=Path, required=True)
    parser.add_argument("--cut-tag-rep1", type=Path, required=True)
    parser.add_argument("--cut-tag-rep2", type=Path, required=True)
    parser.add_argument("--pybigwig-path", type=Path)
    parser.add_argument("--output-table", type=Path, required=True)
    parser.add_argument("--output-sample", type=Path, required=True)
    parser.add_argument("--output-summary", type=Path, required=True)
    return parser.parse_args()


def sample_predictions(path: Path) -> tuple[int, list[dict[str, object]]]:
    rows: list[dict[str, object]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            rows.append(
                {
                    "chrom": row["chrom"],
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                    "candidate_id": row["candidate_id"],
                    "fuzzy_score": float(row["fuzzy_score"]),
                    "hunter_score": float(row["grid_max_raw_score"]),
                    "tss_distance_bp": int(row["tss_distance_bp"]),
                    "tss_context_label": row["tss_context_label"],
                }
            )
    rows.sort(key=lambda row: (float(row["fuzzy_score"]), int(row["start"]), int(row["end"])))
    bins: list[list[dict[str, object]]] = [[] for _ in range(10)]
    for rank, row in enumerate(rows):
        decile = min(9, rank * 10 // len(rows))
        row["fuzzy_score_decile"] = decile + 1
        bins[decile].append(row)
    rng = random.Random(SEED)
    sampled = [row for group in bins for row in rng.sample(group, min(SAMPLE_PER_DECILE, len(group)))]
    sampled.sort(key=lambda row: (int(row["fuzzy_score_decile"]), int(row["start"])))
    return len(rows), sampled


def summarize_group(rows: Sequence[dict[str, object]], decile: int) -> dict[str, object]:
    values = lambda key: [float(row[key]) for row in rows]
    return {
        "fuzzy_score_decile": decile,
        "n": len(rows),
        "fuzzy_score_mean": statistics.fmean(values("fuzzy_score")),
        "fuzzy_score_median": statistics.median(values("fuzzy_score")),
        "zip_ip_mean": statistics.fmean(values("ip_mean")),
        "zip_ip_median": statistics.median(values("ip_mean")),
        "input_mean": statistics.fmean(values("input")),
        "igg_mean": statistics.fmean(values("igg")),
        "log2_zip_background_mean": statistics.fmean(values("log2_zip_background")),
        "log2_zip_background_median": statistics.median(values("log2_zip_background")),
        "cut_tag_mean": statistics.fmean(values("cut_tag_mean")),
        "cut_tag_median": statistics.median(values("cut_tag_mean")),
        "experimental_positive_n": sum(int(row["experimental_positive"]) for row in rows),
        "experimental_positive_rate": rate(rows),
    }


def main() -> None:
    args = parse_args()
    if args.pybigwig_path:
        sys.path.insert(0, str(args.pybigwig_path))
    import pyBigWig  # type: ignore[import-not-found]

    total_candidates, sampled = sample_predictions(args.predictions)
    paths = [
        args.input,
        args.igg,
        args.ip_rep1,
        args.ip_rep2,
        args.cut_tag_rep1,
        args.cut_tag_rep2,
    ]
    tracks = [pyBigWig.open(str(path)) for path in paths]
    try:
        chrom_size = tracks[0].chroms()["Chr1"]
        if any(track.chroms().get("Chr1") != chrom_size for track in tracks):
            raise ValueError("BigWig chromosome sizes do not agree for Chr1")
        for index, row in enumerate(sampled, 1):
            input_signal, igg, ip1, ip2, cut1, cut2 = [
                mean_signal(track, str(row["chrom"]), int(row["start"]), int(row["end"]), chrom_size)
                for track in tracks
            ]
            background = max(input_signal, igg)
            row.update(
                {
                    "input": input_signal,
                    "igg": igg,
                    "ip_rep1": ip1,
                    "ip_rep2": ip2,
                    "ip_mean": (ip1 + ip2) / 2.0,
                    "cut_tag_rep1": cut1,
                    "cut_tag_rep2": cut2,
                    "cut_tag_mean": (cut1 + cut2) / 2.0,
                    "log2_enrichment_rep1": math.log2((ip1 + PSEUDOCOUNT) / (background + PSEUDOCOUNT)),
                    "log2_enrichment_rep2": math.log2((ip2 + PSEUDOCOUNT) / (background + PSEUDOCOUNT)),
                }
            )
            row["log2_zip_background"] = (
                float(row["log2_enrichment_rep1"]) + float(row["log2_enrichment_rep2"])
            ) / 2.0
            if index % 5000 == 0:
                print(f"Measured BigWig signal for {index:,}/{len(sampled):,} sampled candidates")
    finally:
        for track in tracks:
            track.close()

    high_signal_threshold = percentile([float(row["ip_mean"]) for row in sampled], 0.90)
    for row in sampled:
        row["experimental_positive"] = int(
            float(row["ip_rep1"]) > float(row["input"])
            and float(row["ip_rep1"]) > float(row["igg"])
            and float(row["ip_rep2"]) > float(row["input"])
            and float(row["ip_rep2"]) > float(row["igg"])
            and float(row["ip_mean"]) >= high_signal_threshold
        )

    decile_rows = [
        summarize_group([row for row in sampled if row["fuzzy_score_decile"] == decile], decile)
        for decile in range(1, 11)
    ]
    top = [row for row in sampled if row["fuzzy_score_decile"] == 10]
    bottom = [row for row in sampled if row["fuzzy_score_decile"] == 1]
    distal_top = [row for row in top if int(row["tss_distance_bp"]) > 20000]
    distal_bottom = [row for row in bottom if int(row["tss_distance_bp"]) > 20000]
    ci_low, ci_high = bootstrap_rate_difference(top, bottom, random.Random(SEED + 1))
    top_rate, bottom_rate = rate(top), rate(bottom)
    zip_enrichment = [float(row["log2_zip_background"]) for row in sampled]
    cut_tag_signal = [float(row["cut_tag_mean"]) for row in sampled]
    fuzzy_scores = [float(row["fuzzy_score"]) for row in sampled]
    hunter_scores = [float(row["hunter_score"]) for row in sampled]
    summary = {
        "dataset": "GSE252954 (doi:10.1111/pbi.14585)",
        "condition": "Nipponbare rice leaf ZIP-seq; two Z22-IP replicates, input and IgG controls",
        "assembly": "MSU/UGA Rice Genome Annotation Release 7",
        "scope": "chromosome 1",
        "model_status": "frozen Shin-publication preset; no rice tuning",
        "window_bp": WINDOW_BP,
        "sampling_seed": SEED,
        "total_candidates": total_candidates,
        "sampled_candidates": len(sampled),
        "sample_per_decile": SAMPLE_PER_DECILE,
        "pseudocount": PSEUDOCOUNT,
        "experimental_positive_definition": "both ZIP-seq IP replicates > input and IgG, and mean IP in top 10% of sampled windows",
        "high_signal_threshold": high_signal_threshold,
        "spearman_fuzzy_vs_log2_zip_background": spearman(fuzzy_scores, zip_enrichment),
        "spearman_hunter_vs_log2_zip_background": spearman(hunter_scores, zip_enrichment),
        "spearman_fuzzy_vs_cut_tag_mean": spearman(fuzzy_scores, cut_tag_signal),
        "spearman_hunter_vs_cut_tag_mean": spearman(hunter_scores, cut_tag_signal),
        "zip_ip_replicate_pearson": pearson(
            [float(row["ip_rep1"]) for row in sampled], [float(row["ip_rep2"]) for row in sampled]
        ),
        "cut_tag_replicate_pearson": pearson(
            [float(row["cut_tag_rep1"]) for row in sampled],
            [float(row["cut_tag_rep2"]) for row in sampled],
        ),
        "top_decile_positive_rate": top_rate,
        "bottom_decile_positive_rate": bottom_rate,
        "top_vs_bottom_risk_ratio": top_rate / bottom_rate if bottom_rate else math.inf,
        "top_minus_bottom_rate_difference": top_rate - bottom_rate,
        "top_minus_bottom_rate_difference_bootstrap_95ci": [ci_low, ci_high],
        "distal_gt20kb_top_n": len(distal_top),
        "distal_gt20kb_bottom_n": len(distal_bottom),
        "distal_gt20kb_top_positive_rate": rate(distal_top),
        "distal_gt20kb_bottom_positive_rate": rate(distal_bottom),
        "source_files": [path.name for path in paths],
    }

    write_csv(args.output_table, decile_rows, list(decile_rows[0]))
    write_csv(args.output_sample, sampled, list(sampled[0]), gzip_output=True)
    args.output_summary.parent.mkdir(parents=True, exist_ok=True)
    args.output_summary.write_text(json.dumps(summary, indent=2, allow_nan=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, allow_nan=True))


if __name__ == "__main__":
    main()
