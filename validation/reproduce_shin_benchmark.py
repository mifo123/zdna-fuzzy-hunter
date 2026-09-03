#!/usr/bin/env python3
"""Reproduce the reported Shin analysis with the released source and context data.

The script first rebuilds the Z-DNA Hunter part of the analysis-ready feature
table from the released HG Shin FASTA files, audits the 391-to-385 locus
accounting, invokes the public CLI for each operating mode, and verifies the
results against the machine-readable expected values.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path
from typing import Dict


REPO_ROOT = Path(__file__).resolve().parents[1]
DETECTOR = REPO_ROOT / "zdna_fuzzy_detector.py"
BUILDER = REPO_ROOT / "validation" / "build_shin_analysis_features.py"
DATA_DIR = REPO_ROOT / "validation" / "data"
DEFAULT_INPUT = DATA_DIR / "shin_analysis_features.csv"
DEFAULT_CONTEXT = DATA_DIR / "shin_context_features.csv"
SOURCE_DIR = DATA_DIR / "source"
EXPECTED_PATH = DATA_DIR / "shin_expected_metrics.json"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "results" / "shin_reproduction"
MODES = ("balanced", "moderate", "strict")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        help="Use a custom analysis-ready table instead of rebuilding the released one.",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def labelled_metrics(rows: list[dict[str, str]]) -> Dict[str, float | int]:
    counts = {"tp": 0, "fn": 0, "fp": 0, "tn": 0}
    for row in rows:
        label = row.get("shin_label", "").strip().lower()
        if label not in {"positive", "negative"}:
            raise ValueError(f"Unexpected Shin label: {label!r}")
        predicted = int(row["prediction"]) == 1
        if label == "positive" and predicted:
            counts["tp"] += 1
        elif label == "positive":
            counts["fn"] += 1
        elif predicted:
            counts["fp"] += 1
        else:
            counts["tn"] += 1

    recall = counts["tp"] / (counts["tp"] + counts["fn"])
    specificity = counts["tn"] / (counts["tn"] + counts["fp"])
    return {
        **counts,
        "selected": counts["tp"] + counts["fp"],
        "recall": round(recall, 6),
        "specificity": round(specificity, 6),
        "balanced_accuracy": round((recall + specificity) / 2.0, 6),
    }


def verify_metrics(mode: str, actual: Dict[str, float | int], expected: dict[str, object]) -> None:
    expected_mode = expected["modes"][mode]
    mismatches = {
        key: (expected_mode[key], actual.get(key))
        for key in expected_mode
        if actual.get(key) != expected_mode[key]
    }
    if mismatches:
        details = ", ".join(
            f"{key}: expected {wanted}, observed {observed}"
            for key, (wanted, observed) in mismatches.items()
        )
        raise RuntimeError(f"Shin reproduction failed for {mode}: {details}")


def write_summary_csv(rows: list[dict[str, object]], path: Path) -> None:
    columns = [
        "mode",
        "tp",
        "fn",
        "fp",
        "tn",
        "selected",
        "recall",
        "specificity",
        "balanced_accuracy",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    expected = json.loads(EXPECTED_PATH.read_text(encoding="utf-8"))
    output_dir.mkdir(parents=True, exist_ok=True)

    rebuild_performed = args.input is None
    if rebuild_performed:
        if sha256(DEFAULT_CONTEXT) != expected["context_input_sha256"]:
            raise SystemExit("Released author-derived context table checksum is not the expected value")
        for filename, wanted_sha256 in expected["source_sha256"].items():
            if sha256(SOURCE_DIR / filename) != wanted_sha256:
                raise SystemExit(f"Released source FASTA checksum mismatch: {filename}")
        input_path = output_dir / "shin_rebuilt_analysis_features.csv"
        subprocess.run(
            [sys.executable, str(BUILDER), "--output", str(input_path)],
            cwd=REPO_ROOT,
            check=True,
        )
        released_sha256 = sha256(DEFAULT_INPUT)
        rebuilt_sha256 = sha256(input_path)
        if released_sha256 != expected["analysis_input_sha256"]:
            raise SystemExit("Released analysis feature table checksum is not the expected value")
        if rebuilt_sha256 != released_sha256:
            raise SystemExit("Rebuilt analysis feature table differs from the released table")
    else:
        input_path = args.input.resolve()
        if not input_path.exists():
            raise SystemExit(f"Shin feature table not found: {input_path}")

    observed_sha256 = sha256(input_path)
    summary_rows: list[dict[str, object]] = []
    for mode in MODES:
        output_path = output_dir / f"shin_{mode}_predictions.csv"
        command = [
            sys.executable,
            str(DETECTOR),
            "--score-feature-table",
            str(input_path),
            "--preset",
            "shin-publication",
            "--mode",
            mode,
            "--include-all",
            "--output",
            str(output_path),
        ]
        subprocess.run(command, cwd=REPO_ROOT, check=True)
        metrics = labelled_metrics(read_rows(output_path))
        verify_metrics(mode, metrics, expected)
        summary_rows.append({"mode": mode, **metrics})

    summary_csv = output_dir / "shin_benchmark_metrics.csv"
    write_summary_csv(summary_rows, summary_csv)
    summary_json = output_dir / "shin_benchmark_metrics.json"
    summary_json.write_text(
        json.dumps(
            {
                "model_version": expected["model_version"],
                "preset": "shin-publication",
                "rebuild_performed": rebuild_performed,
                "input_kind": expected["input_kind"],
                "input": str(input_path),
                "input_sha256": observed_sha256,
                "loci": expected["loci"],
                "source_loci": expected["source_loci"],
                "excluded_loci": expected["excluded_loci"],
                "modes": {row["mode"]: {k: v for k, v in row.items() if k != "mode"} for row in summary_rows},
                "matches_released_expected_values": True,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"Verified all {len(MODES)} Shin operating modes after source/cohort audit.")
    print(f"Metrics: {summary_csv}")
    print(f"Machine-readable summary: {summary_json}")


if __name__ == "__main__":
    main()
