#!/usr/bin/env python3
"""Independent U2OS Z-DNA ChIP-seq check for fixed ZDNA-Fuzzy Hunter calls.

The analysis uses deterministic stratified sampling across fuzzy-score deciles.
For each sampled candidate it measures exact mean BigWig coverage in a fixed
500-bp window centred on the prediction.  A conservative experimental-positive
label requires (i) Z22-IP coverage greater than the matched input in both WT
replicates and (ii) mean Z22-IP coverage in the top 10% of sampled windows.
No U2OS value is used to tune the model, its threshold, or its fuzzy rules.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import random
import statistics
import sys
from pathlib import Path
from typing import Iterable, Sequence


SEED = 20260803
SAMPLE_PER_DECILE = 2000
WINDOW_BP = 500
PSEUDOCOUNT = 1.0
BOOTSTRAPS = 2000


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--input-rep1", type=Path, required=True)
    parser.add_argument("--input-rep2", type=Path, required=True)
    parser.add_argument("--ip-rep1", type=Path, required=True)
    parser.add_argument("--ip-rep2", type=Path, required=True)
    parser.add_argument("--pybigwig-path", type=Path)
    parser.add_argument("--output-table", type=Path, required=True)
    parser.add_argument("--output-sample", type=Path, required=True)
    parser.add_argument("--output-summary", type=Path, required=True)
    return parser.parse_args()


def percentile(values: Sequence[float], q: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return math.nan
    position = (len(ordered) - 1) * q
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def average_ranks(values: Sequence[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    index = 0
    while index < len(order):
        end = index + 1
        while end < len(order) and values[order[end]] == values[order[index]]:
            end += 1
        rank = (index + 1 + end) / 2.0
        for ordered_index in order[index:end]:
            ranks[ordered_index] = rank
        index = end
    return ranks


def pearson(left: Sequence[float], right: Sequence[float]) -> float:
    left_mean = statistics.fmean(left)
    right_mean = statistics.fmean(right)
    numerator = sum((x - left_mean) * (y - right_mean) for x, y in zip(left, right))
    left_ss = sum((x - left_mean) ** 2 for x in left)
    right_ss = sum((y - right_mean) ** 2 for y in right)
    denominator = math.sqrt(left_ss * right_ss)
    return numerator / denominator if denominator else math.nan


def spearman(left: Sequence[float], right: Sequence[float]) -> float:
    return pearson(average_ranks(left), average_ranks(right))


def rate(rows: Sequence[dict[str, object]]) -> float:
    return statistics.fmean(float(row["experimental_positive"]) for row in rows) if rows else math.nan


def bootstrap_rate_difference(
    top: Sequence[dict[str, object]], bottom: Sequence[dict[str, object]], rng: random.Random
) -> tuple[float, float]:
    differences: list[float] = []
    for _ in range(BOOTSTRAPS):
        top_rate = statistics.fmean(
            float(top[rng.randrange(len(top))]["experimental_positive"]) for _ in range(len(top))
        )
        bottom_rate = statistics.fmean(
            float(bottom[rng.randrange(len(bottom))]["experimental_positive"]) for _ in range(len(bottom))
        )
        differences.append(top_rate - bottom_rate)
    return percentile(differences, 0.025), percentile(differences, 0.975)


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


def mean_signal(bigwig: object, chrom: str, start: int, end: int, chrom_size: int) -> float:
    midpoint = (start + end) // 2
    half = WINDOW_BP // 2
    window_start = max(0, midpoint - half)
    window_end = min(chrom_size, midpoint + half)
    value = bigwig.stats(chrom, window_start, window_end, type="mean", exact=True)[0]
    return 0.0 if value is None or not math.isfinite(value) else float(value)


def summarize_group(rows: Sequence[dict[str, object]], decile: int) -> dict[str, object]:
    values = lambda key: [float(row[key]) for row in rows]
    return {
        "fuzzy_score_decile": decile,
        "n": len(rows),
        "fuzzy_score_mean": statistics.fmean(values("fuzzy_score")),
        "fuzzy_score_median": statistics.median(values("fuzzy_score")),
        "z22_ip_mean": statistics.fmean(values("ip_mean")),
        "z22_ip_median": statistics.median(values("ip_mean")),
        "input_mean": statistics.fmean(values("input_mean")),
        "input_median": statistics.median(values("input_mean")),
        "log2_ip_input_mean": statistics.fmean(values("log2_ip_input")),
        "log2_ip_input_median": statistics.median(values("log2_ip_input")),
        "experimental_positive_n": sum(int(row["experimental_positive"]) for row in rows),
        "experimental_positive_rate": rate(rows),
    }


def write_csv(path: Path, rows: Iterable[dict[str, object]], fieldnames: Sequence[str], gzip_output: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if gzip_output else open
    with opener(path, "wt", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    if args.pybigwig_path:
        sys.path.insert(0, str(args.pybigwig_path))
    import pyBigWig  # type: ignore[import-not-found]

    total_candidates, sampled = sample_predictions(args.predictions)
    paths = [args.input_rep1, args.input_rep2, args.ip_rep1, args.ip_rep2]
    tracks = [pyBigWig.open(str(path)) for path in paths]
    try:
        chrom_size = tracks[0].chroms()["chr22"]
        if any(track.chroms().get("chr22") != chrom_size for track in tracks):
            raise ValueError("BigWig chromosome sizes do not agree for chr22")
        for index, row in enumerate(sampled, 1):
            signals = [
                mean_signal(track, str(row["chrom"]), int(row["start"]), int(row["end"]), chrom_size)
                for track in tracks
            ]
            input1, input2, ip1, ip2 = signals
            row.update(
                {
                    "input_rep1": input1,
                    "input_rep2": input2,
                    "ip_rep1": ip1,
                    "ip_rep2": ip2,
                    "input_mean": (input1 + input2) / 2.0,
                    "ip_mean": (ip1 + ip2) / 2.0,
                    "log2_enrichment_rep1": math.log2((ip1 + PSEUDOCOUNT) / (input1 + PSEUDOCOUNT)),
                    "log2_enrichment_rep2": math.log2((ip2 + PSEUDOCOUNT) / (input2 + PSEUDOCOUNT)),
                }
            )
            row["log2_ip_input"] = (
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
            float(row["ip_rep1"]) > float(row["input_rep1"])
            and float(row["ip_rep2"]) > float(row["input_rep2"])
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
    rng = random.Random(SEED + 1)
    ci_low, ci_high = bootstrap_rate_difference(top, bottom, rng)
    top_rate, bottom_rate = rate(top), rate(bottom)
    risk_ratio = top_rate / bottom_rate if bottom_rate else math.inf
    fuzzy_scores = [float(row["fuzzy_score"]) for row in sampled]
    enrichments = [float(row["log2_ip_input"]) for row in sampled]
    summary = {
        "dataset": "GSE290662",
        "condition": "U2OS WT (raw_CT), Z22-IP versus matched input, two biological replicates",
        "assembly": "hg19",
        "scope": "chromosome 22",
        "model_status": "frozen Shin-publication preset; no U2OS tuning",
        "window_bp": WINDOW_BP,
        "sampling_seed": SEED,
        "total_candidates": total_candidates,
        "sampled_candidates": len(sampled),
        "sample_per_decile": SAMPLE_PER_DECILE,
        "pseudocount": PSEUDOCOUNT,
        "experimental_positive_definition": "IP > matched input in both replicates and mean IP in top 10% of sampled windows",
        "high_signal_threshold": high_signal_threshold,
        "spearman_fuzzy_vs_log2_ip_input": spearman(fuzzy_scores, enrichments),
        "spearman_hunter_vs_log2_ip_input": spearman(
            [float(row["hunter_score"]) for row in sampled], enrichments
        ),
        "replicate_log2_enrichment_pearson": pearson(
            [float(row["log2_enrichment_rep1"]) for row in sampled],
            [float(row["log2_enrichment_rep2"]) for row in sampled],
        ),
        "top_decile_positive_rate": top_rate,
        "bottom_decile_positive_rate": bottom_rate,
        "top_vs_bottom_risk_ratio": risk_ratio,
        "top_minus_bottom_rate_difference": top_rate - bottom_rate,
        "top_minus_bottom_rate_difference_bootstrap_95ci": [ci_low, ci_high],
        "distal_gt20kb_top_n": len(distal_top),
        "distal_gt20kb_bottom_n": len(distal_bottom),
        "distal_gt20kb_top_positive_rate": rate(distal_top),
        "distal_gt20kb_bottom_positive_rate": rate(distal_bottom),
        "source_files": [path.name for path in paths],
    }

    table_fields = list(decile_rows[0])
    write_csv(args.output_table, decile_rows, table_fields)
    sample_fields = list(sampled[0])
    write_csv(args.output_sample, sampled, sample_fields, gzip_output=True)
    args.output_summary.parent.mkdir(parents=True, exist_ok=True)
    args.output_summary.write_text(json.dumps(summary, indent=2, allow_nan=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, allow_nan=True))


if __name__ == "__main__":
    main()
