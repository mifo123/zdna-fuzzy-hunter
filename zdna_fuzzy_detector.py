#!/usr/bin/env python3
"""Genome-scale Z-DNA candidate prioritization with Z-DNA Hunter and fuzzy logic.

The default backend runs the Z-DNA Hunter core locally.  An optional IBP API
backend is retained for compatibility with hosted analyses.  Both backends
produce 0-based, half-open intervals that are merged, annotated with nearest-TSS
context and scored by the same fixed fuzzy expert model.
"""

from __future__ import annotations

import argparse
import csv
import getpass
import json
import math
import os
import platform
import re
import statistics
import time
from bisect import bisect_left
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


Number = float
MODEL_VERSION = "1.1.0"


def hunter_config(config_id: str, min_sequence_size: int, threshold: Number) -> Dict[str, object]:
    return {
        "config_id": config_id,
        "model": "model2",
        "min_sequence_size": min_sequence_size,
        "threshold": threshold,
        "score_gc": 2.0,
        "score_gtac": 1.0,
        "score_at": 0.5,
        "score_oth": 0.0,
    }


HUNTER_PRESETS = {
    "genome-balanced": [
        hunter_config("m2_l8_t30", 8, 30.0),
        hunter_config("m2_l10_t60", 10, 60.0),
    ],
    "genome-strict": [
        hunter_config("m2_l10_t60", 10, 60.0),
    ],
    "shin-publication": [
        hunter_config("m2_l6_t30", 6, 30.0),
        hunter_config("m2_l10_t60", 10, 60.0),
    ],
}

HUNTER_CONFIGS = list(HUNTER_PRESETS["genome-balanced"])

FUZZY_THRESHOLD = 19.157
DEFAULT_GENOME_MIN_SCORE = 30.0
MODERATE_TSS_BP = 10
MODERATE_MAX_SCORE = 19.5
STRICT_TSS_BP = 50
STRICT_MAX_SIGNAL = 38.51
STRICT_MAX_OVERLAP_PCT = 2.0


@dataclass(frozen=True)
class FastaRecord:
    name: str
    sequence: str
    output_chrom: str
    coordinate_offset: int = 0


@dataclass(frozen=True)
class Interval:
    """A genomic interval using 0-based, half-open ``[start, end)`` coordinates."""

    chrom: str
    start: int
    end: int
    score: Number
    source: str

    @property
    def length(self) -> int:
        return max(0, self.end - self.start)


@dataclass(frozen=True)
class AnnotationTrack:
    kind: str
    name: str
    intervals_by_chrom: Dict[str, List[Tuple[int, int]]]


@dataclass(frozen=True)
class TunedParams:
    """Fixed publication parameters selected by CV-guided Shin random search.

    The same folds contributed to hyperparameter selection and to the reported
    cross-validated estimate.  Those results are therefore tuning-aware internal
    validation, not a fully nested or independent estimate of generalization.
    """

    signal_weighted: Number = 0.37921
    signal_zdna: Number = 0.051894
    signal_vote: Number = 0.198225
    signal_model1: Number = 0.191368
    signal_model2: Number = 0.179302
    signal_all5_bonus: Number = 3.089634
    context_tss: Number = 0.47675
    context_promoter: Number = 0.02146
    context_marks: Number = 0.408899
    evidence_source: Number = 0.459564
    evidence_model1: Number = 0.130052
    evidence_model2: Number = 0.355236
    evidence_balance: Number = 0.055148
    evidence_all5_bonus: Number = 5.240791
    bias_repeat: Number = 0.296978
    bias_disagreement: Number = 0.16273
    bias_model2_only: Number = 0.284105
    bias_no_hit: Number = 0.256187
    final_signal: Number = 0.439786
    final_evidence: Number = 0.448789
    final_context: Number = 0.010134
    final_feasibility: Number = 0.101291
    final_bias: Number = 0.471126
    rule_mix: Number = 0.149236
    offset: Number = 4.36488


class FuzzyTerm:
    def __init__(self, name: str, points: Tuple[Number, Number, Number, Number]) -> None:
        self.name = name
        self.a, self.b, self.c, self.d = points

    def membership(self, value: Number) -> Number:
        if value <= self.a:
            return 1.0 if self.a == self.b and value == self.a else 0.0
        if self.a < value < self.b:
            return (value - self.a) / max(self.b - self.a, 1e-12)
        if self.b <= value <= self.c:
            return 1.0
        if self.c < value < self.d:
            return (self.d - value) / max(self.d - self.c, 1e-12)
        if value >= self.d:
            return 1.0 if self.c == self.d and value == self.d else 0.0
        return 0.0


TERMS = {
    "very_low": FuzzyTerm("very_low", (0, 0, 8, 22)),
    "low": FuzzyTerm("low", (12, 24, 36, 50)),
    "medium": FuzzyTerm("medium", (38, 50, 62, 74)),
    "high": FuzzyTerm("high", (62, 74, 84, 94)),
    "very_high": FuzzyTerm("very_high", (86, 94, 100, 100)),
}

BIAS_TERMS = {
    "very_low": FuzzyTerm("very_low", (0, 0, 6, 16)),
    "low": FuzzyTerm("low", (8, 18, 28, 40)),
    "medium": FuzzyTerm("medium", (30, 42, 55, 68)),
    "high": FuzzyTerm("high", (58, 70, 82, 92)),
    "very_high": FuzzyTerm("very_high", (84, 94, 100, 100)),
}

FUZZY_RULE_SPECS = [
    {"id": "AF1", "if": "signal is very_high AND evidence is very_high AND bias is very_low", "then": 98.0},
    {"id": "AF2", "if": "signal is very_high AND evidence is high_or_very_high AND bias is very_low_or_low", "then": 93.0},
    {"id": "AF3", "if": "signal is high AND evidence is very_high AND context is medium_or_higher AND bias is very_low_or_low", "then": 89.0},
    {"id": "AF4", "if": "signal is high AND evidence is high AND bias is low_or_medium", "then": 82.0},
    {"id": "AF5", "if": "signal is medium AND evidence is high_or_very_high AND context is high", "then": 72.0},
    {"id": "AF6", "if": "signal is high AND evidence is medium AND context is low", "then": 64.0},
    {"id": "AF7", "if": "signal is medium AND evidence is medium AND bias is low", "then": 58.0},
    {"id": "AF8", "if": "signal is very_low AND evidence is low", "then": 12.0},
    {"id": "AF9", "if": "bias is high_or_very_high AND evidence is low_or_medium", "then": 22.0},
    {"id": "AF10", "if": "signal is high AND cross-model balance is low", "then": 48.0},
    {"id": "AF11", "if": "signal is very_high AND cross-model balance is very_high AND feasibility is medium", "then": 91.0},
    {"id": "AF12", "if": "signal is low AND bias is high_or_very_high", "then": 14.0},
]

MODEL2_RUN_PATTERN = re.compile(r"[AG][CT](?:[AG][CT])*[AG]?|[CT][AG](?:[CT][AG])*[CT]?")


@dataclass(frozen=True)
class LocalHunterHit:
    """One local Z-DNA Hunter hit in sequence-local half-open coordinates."""

    start: int
    end: int
    score_percent: Number


def hunter_pair_score(first: str, second: str, config: Dict[str, object]) -> Number:
    """Return the configured Z-DNA Hunter score for one adjacent base pair."""
    pair = (first.upper(), second.upper())
    if pair in {("G", "C"), ("C", "G")}:
        return parse_float(config["score_gc"])
    if pair in {("G", "T"), ("T", "G"), ("A", "C"), ("C", "A")}:
        return parse_float(config["score_gtac"])
    if pair in {("A", "T"), ("T", "A")}:
        return parse_float(config["score_at"])
    return parse_float(config.get("score_oth"))


def positive_score_runs(sequence: str, config: Dict[str, object]) -> Iterable[Tuple[int, str, Number]]:
    """Yield maximal runs whose adjacent pairs have a positive Hunter score."""
    sequence = sequence.upper()
    is_standard_model2 = (
        str(config.get("model", "")).lower() == "model2"
        and parse_float(config.get("score_gc")) > 0.0
        and parse_float(config.get("score_gtac")) > 0.0
        and parse_float(config.get("score_at")) > 0.0
        and parse_float(config.get("score_oth")) == 0.0
    )
    if is_standard_model2:
        for match in MODEL2_RUN_PATTERN.finditer(sequence):
            run = match.group(0)
            score_sum = sum(hunter_pair_score(run[index], run[index + 1], config) for index in range(len(run) - 1))
            yield match.start(), run, score_sum
        return

    run_length = 1
    score_sum = 0.0
    for index in range(len(sequence)):
        pair_score = hunter_pair_score(sequence[index], sequence[index + 1], config) if index < len(sequence) - 1 else 0.0
        if pair_score > 0.0:
            run_length += 1
            score_sum += pair_score
            continue
        if run_length > 1:
            start = index - run_length + 1
            yield start, sequence[start : start + run_length], score_sum
        run_length = 1
        score_sum = 0.0


def local_hunter_hits(sequence: str, config: Dict[str, object]) -> Iterable[LocalHunterHit]:
    """Run the local Java-compatible Z-DNA Hunter state machine on a sequence."""
    minimum_length = int(config["min_sequence_size"])
    threshold = parse_float(config["threshold"])
    max_pair_score = max(
        parse_float(config["score_gc"]),
        parse_float(config["score_gtac"]),
        parse_float(config["score_at"]),
    )
    for start, run, score_sum in positive_score_runs(sequence, config):
        if len(run) < minimum_length:
            continue
        score = score_sum / 2.0
        max_possible = ((len(run) - 1) * max_pair_score) / 2.0
        score_percent = 100.0 * score / max_possible if max_possible > 0.0 else 0.0
        if score_percent >= threshold:
            yield LocalHunterHit(start=start, end=start + len(run), score_percent=score_percent)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run local or API-backed Z-DNA Hunter and fuzzy prioritization on a FASTA genome, chromosome or region."
    )
    parser.add_argument("--fasta", type=Path, default=None, help="Input FASTA file.")
    parser.add_argument("--tss", type=Path, default=None, help="TSS annotation in SGA, BED or CSV/TSV format.")
    parser.add_argument("--output", type=Path, default=None, help="Output CSV or bedGraph path.")
    parser.add_argument("--format", choices=["csv", "bedgraph"], default=None, help="Output format; inferred by extension.")
    parser.add_argument("--mode", choices=["balanced", "moderate", "strict"], default="balanced")
    parser.add_argument(
        "--preset",
        choices=sorted(HUNTER_PRESETS),
        default="genome-balanced",
        help="Z-DNA Hunter preset. Default is practical genome-scale filtering.",
    )
    parser.add_argument(
        "--min-score",
        type=float,
        default=None,
        help="Minimum fuzzy score for final calls. Default: 30 for genome presets, publication threshold for shin-publication.",
    )
    parser.add_argument(
        "--hunter-config",
        type=Path,
        default=None,
        help="Optional CSV with custom Z-DNA Hunter configurations. Overrides --preset configs for this run.",
    )
    parser.add_argument(
        "--hunter-backend",
        choices=["local", "api"],
        default="local",
        help="Run the bundled local Hunter core (default) or submit scans to the IBP API.",
    )
    parser.add_argument("--email", default=os.getenv("IBP_EMAIL"), help="IBP account email; used only with --hunter-backend api.")
    parser.add_argument("--password", default=os.getenv("IBP_PASSWORD"), help="IBP account password; used only with --hunter-backend api.")
    parser.add_argument("--server", default=None, help="Optional IBP API server URL.")
    parser.add_argument("--work-dir", type=Path, default=Path("zdna_fuzzy_runs"), help="Intermediate files directory.")
    parser.add_argument("--run-name", default=None, help="Reusable API tag/name prefix. Default is based on FASTA stem.")
    parser.add_argument("--chromosomes", default="", help="Comma-separated FASTA records to process; default is all.")
    parser.add_argument("--region", default="", help="Optional 1-based closed interval, for example chr2L:1000-5000.")
    parser.add_argument("--coordinate-offset", type=int, default=0, help="Offset added to FASTA-local coordinates.")
    parser.add_argument(
        "--api-coordinate-base",
        choices=["0", "1"],
        default="0",
        help="Coordinate base used in Z-DNA Hunter exports. Default follows previous IBP exports.",
    )
    parser.add_argument(
        "--tss-coordinate-base",
        choices=["auto", "0", "1"],
        default="auto",
        help="Coordinate base of point TSS annotations. Auto treats SGA as 1-based and BED as 0-based.",
    )
    parser.add_argument("--promoter-window", type=int, default=1000, help="Nearest-TSS context window in bp.")
    parser.add_argument(
        "--annotation-bed",
        action="append",
        default=[],
        metavar="NAME=PATH",
        help=(
            "Optional genomic-context BED track for candidate overlap annotation, "
            "for example promoters=tracks/promoters.bed. Can be passed multiple times."
        ),
    )
    parser.add_argument(
        "--epigenomic-bed",
        action="append",
        default=[],
        metavar="NAME=PATH",
        help=(
            "Optional BED/narrowPeak/broadPeak track for candidate overlap annotation, "
            "for example ATAC=tracks/atac.bed. Can be passed multiple times."
        ),
    )
    parser.add_argument("--poll-seconds", type=float, default=2.0, help="Polling interval for API jobs.")
    parser.add_argument("--reuse", action=argparse.BooleanOptionalAction, default=True, help="Reuse API objects with matching tags.")
    parser.add_argument("--keep-intermediate", action="store_true", help="Keep uploaded FASTA slices and raw exports.")
    parser.add_argument("--include-all", action="store_true", help="Write all candidates, not only selected positives.")
    parser.add_argument("--summary-output", type=Path, default=None, help="Optional JSON summary path.")
    validation = parser.add_argument_group("validation / reproducibility")
    validation.add_argument(
        "--require-grid-hit",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="In feature-table mode, gate candidates without Z-DNA Hunter grid support.",
    )
    validation.add_argument(
        "--score-feature-table",
        dest="feature_table",
        type=Path,
        default=None,
        help="Validation mode: score an already prepared candidate feature table and skip API execution.",
    )
    validation.add_argument("--feature-table", dest="feature_table", type=Path, help=argparse.SUPPRESS)
    inspection = parser.add_argument_group("model inspection")
    inspection.add_argument(
        "--describe-model",
        action="store_true",
        help="Print the complete fixed model specification as JSON and exit.",
    )
    inspection.add_argument(
        "--model-spec-output",
        type=Path,
        default=None,
        help="Optional path receiving the complete model specification JSON.",
    )
    return parser.parse_args()


def clamp(value: Number, lower: Number = 0.0, upper: Number = 100.0) -> Number:
    return max(lower, min(upper, value))


def load_hunter_config(path: Path) -> List[Dict[str, object]]:
    required = {
        "config_id",
        "model",
        "min_sequence_size",
        "threshold",
        "score_gc",
        "score_gtac",
        "score_at",
        "score_oth",
    }
    rows = read_delimited_text(path.read_text(encoding="utf-8-sig", errors="replace"))
    if not rows:
        raise SystemExit(f"Hunter config file is empty: {path}")
    missing = required - set(rows[0])
    if missing:
        raise SystemExit(f"Hunter config file is missing columns: {', '.join(sorted(missing))}")
    configs: List[Dict[str, object]] = []
    seen_ids: set[str] = set()
    for row in rows:
        config_id = str(row["config_id"]).strip()
        model = str(row["model"]).strip()
        if not config_id:
            raise SystemExit("Hunter config contains an empty config_id.")
        if config_id in seen_ids:
            raise SystemExit(f"Hunter config contains duplicate config_id: {config_id}")
        seen_ids.add(config_id)
        if model not in {"model1", "model2"}:
            raise SystemExit(f"Unsupported Hunter model for {config_id}: {model}")
        min_sequence_size = parse_int(row.get("min_sequence_size"))
        threshold = parse_float(row.get("threshold"))
        if min_sequence_size < 1:
            raise SystemExit(f"Invalid min_sequence_size for {config_id}: {min_sequence_size}")
        if threshold <= 0.0 or threshold > 100.0:
            raise SystemExit(f"Invalid threshold for {config_id}: {threshold}")
        configs.append(
            {
                "config_id": config_id,
                "model": model,
                "min_sequence_size": min_sequence_size,
                "threshold": threshold,
                "score_gc": parse_float(row.get("score_gc")),
                "score_gtac": parse_float(row.get("score_gtac")),
                "score_at": parse_float(row.get("score_at")),
                "score_oth": parse_float(row.get("score_oth")),
            }
        )
    return configs


def apply_preset(args: argparse.Namespace) -> None:
    global HUNTER_CONFIGS
    if args.hunter_config is not None:
        HUNTER_CONFIGS = load_hunter_config(args.hunter_config)
    else:
        HUNTER_CONFIGS = [dict(config) for config in HUNTER_PRESETS[args.preset]]


def effective_min_score(args: argparse.Namespace) -> Number:
    if args.min_score is not None:
        return float(args.min_score)
    if args.feature_table is not None or args.preset == "shin-publication":
        return FUZZY_THRESHOLD
    return DEFAULT_GENOME_MIN_SCORE


def slug(value: str) -> str:
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    return text.strip("_") or "sequence"


def canonical_chrom(value: object) -> str:
    text = str(value).strip()
    return text[3:].lower() if text.lower().startswith("chr") else text.lower()


def split_csv_arg(value: str) -> List[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def parse_float(value: object, default: Number = 0.0) -> Number:
    if value is None:
        return default
    text = str(value).strip()
    if not text:
        return default
    text = text.replace("%", "").replace(",", ".")
    try:
        return float(text)
    except ValueError:
        return default


def parse_int(value: object, default: int = 0) -> int:
    return int(round(parse_float(value, float(default))))


def parse_bool(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def detect_column(fieldnames: Sequence[str], candidates: Sequence[str]) -> Optional[str]:
    by_norm = {re.sub(r"[^a-z0-9]+", "", name.lower()): name for name in fieldnames}
    for candidate in candidates:
        norm = re.sub(r"[^a-z0-9]+", "", candidate.lower())
        if norm in by_norm:
            return by_norm[norm]
    return None


def read_delimited_text(text: str) -> List[Dict[str, str]]:
    if not text.strip():
        return []
    sample = text[:4096]
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",;\t")
    except csv.Error:
        dialect = csv.excel
    return list(csv.DictReader(text.splitlines(), dialect=dialect))


def write_csv(rows: Sequence[Dict[str, object]], path: Path, columns: Optional[Sequence[str]] = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if columns is not None:
        fieldnames = list(columns)
    else:
        fieldnames = output_columns()
        for row in rows:
            for key in row:
                if key not in fieldnames:
                    fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def first_present(row: Dict[str, str], names: Sequence[str], default: object = "") -> object:
    for name in names:
        value = row.get(name)
        if value is not None and str(value).strip() != "":
            return value
    return default


def rename_derived_feature_columns(row: Dict[str, object]) -> Dict[str, object]:
    renamed = dict(row)
    aliases = [
        ("promoter_overlap", "tss_proximal"),
        ("regulatory_marks", "context_support_bin"),
        ("repeat_overlap_pct", "heuristic_bias_score"),
        ("primer_uniqueness", "candidate_uniqueness_score"),
    ]
    for old, new in aliases:
        if old in renamed and new not in renamed:
            renamed[new] = renamed[old]
        renamed.pop(old, None)
    return renamed


def read_fasta(path: Path) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    name: Optional[str] = None
    parts: List[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(parts).upper()))
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(re.sub(r"[^A-Za-z]", "", line))
    if name is not None:
        records.append((name, "".join(parts).upper()))
    if not records:
        raise SystemExit(f"No FASTA records found in {path}")
    return records


def parse_region(region: str) -> Optional[Tuple[str, int, int]]:
    if not region:
        return None
    match = re.match(r"^([^:]+):(\d+)-(\d+)$", region.replace(",", ""))
    if not match:
        raise SystemExit("--region must look like chrom:start-end")
    chrom = match.group(1)
    start_1 = int(match.group(2))
    end_1 = int(match.group(3))
    if start_1 < 1 or end_1 < start_1:
        raise SystemExit("--region must use a valid 1-based closed interval")
    return chrom, start_1, end_1


def selected_fasta_records(args: argparse.Namespace) -> List[FastaRecord]:
    raw_records = read_fasta(args.fasta)
    requested = {canonical_chrom(item) for item in split_csv_arg(args.chromosomes)}
    region = parse_region(args.region)
    selected: List[FastaRecord] = []
    for name, sequence in raw_records:
        if requested and canonical_chrom(name) not in requested:
            continue
        output_chrom = name
        offset = args.coordinate_offset
        seq = sequence
        if region:
            region_chrom, start_1, end_1 = region
            if canonical_chrom(name) != canonical_chrom(region_chrom):
                continue
            seq = sequence[start_1 - 1 : end_1]
            output_chrom = region_chrom
            offset += start_1 - 1
        selected.append(FastaRecord(name=name, sequence=seq, output_chrom=output_chrom, coordinate_offset=offset))
    if not selected:
        raise SystemExit("No FASTA records matched --chromosomes/--region.")
    return selected


def run_local_hunter(records: Sequence[FastaRecord]) -> Dict[str, List[Interval]]:
    """Run every configured Hunter scan locally and return half-open intervals."""
    intervals_by_chrom: Dict[str, List[Interval]] = {}
    for record in records:
        chrom_intervals = intervals_by_chrom.setdefault(record.output_chrom, [])
        for config in HUNTER_CONFIGS:
            config_id = str(config["config_id"])
            hit_count = 0
            for hit in local_hunter_hits(record.sequence, config):
                chrom_intervals.append(
                    Interval(
                        chrom=record.output_chrom,
                        start=record.coordinate_offset + hit.start,
                        end=record.coordinate_offset + hit.end,
                        score=hit.score_percent,
                        source=config_id,
                    )
                )
                hit_count += 1
            print(f"Local Hunter {record.output_chrom}/{config_id}: {hit_count} interval(s)", flush=True)
        chrom_intervals.sort(key=lambda item: (item.start, item.end, item.source))
    return intervals_by_chrom


def looks_like_header(parts: Sequence[str]) -> bool:
    if (
        len(parts) >= 4
        and parts[1].strip().lower() == "tss"
        and re.fullmatch(r"-?\d+", parts[2].strip())
        and parts[3].strip() in {"+", "-"}
    ):
        return False
    header_tokens = {
        "chrom",
        "chromosome",
        "chr",
        "seqname",
        "seqid",
        "start",
        "end",
        "strand",
        "position",
        "pos",
        "tss",
    }
    normalized = {re.sub(r"[^a-z0-9]+", "", part.lower()) for part in parts}
    return bool(header_tokens & normalized)


def parse_tss_no_header(parts: Sequence[str], coordinate_base: Optional[int]) -> Optional[Tuple[str, int]]:
    """Parse SGA, BED-like or point TSS records into a 0-based point coordinate."""
    if len(parts) < 2:
        return None
    chrom = parts[0]
    if (
        len(parts) >= 4
        and parts[1].strip().lower() == "tss"
        and re.fullmatch(r"-?\d+", parts[2])
        and parts[3] in {"+", "-"}
    ):
        base = 1 if coordinate_base is None else coordinate_base
        return chrom, int(parts[2]) - base
    if len(parts) >= 3 and re.fullmatch(r"[+-]", parts[2]) and re.fullmatch(r"-?\d+", parts[1]):
        base = 0 if coordinate_base is None else coordinate_base
        pos = int(parts[1]) - base
        return chrom, pos
    if len(parts) >= 3 and re.fullmatch(r"-?\d+", parts[1]) and re.fullmatch(r"-?\d+", parts[2]):
        start = int(parts[1])
        end = int(parts[2])
        strand = next((part for part in parts[3:8] if part in {"+", "-"}), "+")
        pos = max(start, end - 1) if strand == "-" else start
        return chrom, pos
    if re.fullmatch(r"-?\d+", parts[1]):
        base = 0 if coordinate_base is None else coordinate_base
        pos = int(parts[1]) - base
        return chrom, pos
    return None


def load_tss(path: Path, coordinate_base: object = "auto") -> Dict[str, List[int]]:
    """Load TSS annotations, normalizing all supported formats to 0-based points."""
    requested_base = None if str(coordinate_base) == "auto" else int(str(coordinate_base))
    text = path.read_text(encoding="utf-8-sig", errors="replace")
    data: Dict[str, List[int]] = {}
    non_comment = [line for line in text.splitlines() if line.strip() and not line.lstrip().startswith("#")]
    if not non_comment:
        return data
    first_parts = re.split(r"[\t,; ]+", non_comment[0].strip())
    if looks_like_header(first_parts):
        if not any(delimiter in non_comment[0] for delimiter in ("\t", ",", ";")):
            header = first_parts
            rows = [
                dict(zip(header, re.split(r"[\t,; ]+", line.strip())))
                for line in non_comment[1:]
                if line.strip()
            ]
        else:
            rows = read_delimited_text("\n".join(non_comment))
        for row in rows:
            fields = list(row.keys())
            chrom_col = detect_column(fields, ["chrom", "chromosome", "chr", "seqname", "seqid"])
            pos_col = detect_column(fields, ["tss", "position", "pos", "site"])
            start_col = detect_column(fields, ["start", "tx_start", "transcription_start"])
            end_col = detect_column(fields, ["end", "tx_end", "transcription_end"])
            strand_col = detect_column(fields, ["strand", "orientation"])
            if chrom_col is None:
                continue
            chrom = str(row.get(chrom_col, "")).strip()
            if not chrom:
                continue
            if pos_col is not None:
                base = 0 if requested_base is None else requested_base
                pos = parse_int(row.get(pos_col)) - base
            elif start_col is not None:
                start = parse_int(row.get(start_col))
                end = parse_int(row.get(end_col), start) if end_col else start
                strand = str(row.get(strand_col, "+")).strip() if strand_col else "+"
                pos = max(start, end - 1) if strand == "-" else start
            else:
                continue
            data.setdefault(canonical_chrom(chrom), []).append(pos)
    else:
        for line in non_comment:
            parts = re.split(r"[\t,; ]+", line.strip())
            parsed = parse_tss_no_header(parts, requested_base)
            if parsed is None:
                continue
            chrom, pos = parsed
            data.setdefault(canonical_chrom(chrom), []).append(pos)
    for positions in data.values():
        positions.sort()
    return data


def nearest_distance(positions: Sequence[int], start: int, end: int) -> Optional[int]:
    if not positions:
        return None
    mid = (start + end) // 2
    insertion = bisect_left(positions, mid)
    best: Optional[int] = None
    for idx in (insertion - 1, insertion, insertion + 1):
        if 0 <= idx < len(positions):
            pos = positions[idx]
            if start <= pos < end:
                distance = 0
            elif pos < start:
                distance = start - pos
            else:
                distance = pos - end
            best = distance if best is None else min(best, distance)
    return best


def prompt_credentials(args: argparse.Namespace) -> Tuple[str, str]:
    email = args.email or input("IBP email: ").strip()
    password = args.password or getpass.getpass("IBP password: ")
    if not email or not password:
        raise SystemExit("IBP email and password are required.")
    return email, password


def load_ibp_modules():
    try:
        from DNA_analyser_IBP.api import Api
        from DNA_analyser_IBP.adapters.validations import ApiEmptyResponse
        from DNA_analyser_IBP.type import Types
        from DNA_analyser_IBP.utils import normalize_name
    except ImportError as exc:
        raise SystemExit(
            "Missing DNA_analyser_IBP. Install the optional API requirements first:\n"
            "  pip install -r requirements-api.txt"
        ) from exc
    return Api, ApiEmptyResponse, Types, normalize_name


def make_api(args: argparse.Namespace):
    Api, _, _, _ = load_ibp_modules()
    email, password = prompt_credentials(args)
    kwargs = {"email": email, "password": password}
    if args.server:
        kwargs["server"] = args.server
    return Api(**kwargs)


def ports(api):
    return getattr(api, "_Api__ports")


def object_tags(row: Dict[str, object]) -> List[str]:
    return [item.strip() for item in str(row.get("tags", "")).split(",") if item.strip()]


def sequence_records(api, tags: Sequence[str]) -> List[Dict[str, object]]:
    _, ApiEmptyResponse, _, _ = load_ibp_modules()
    try:
        return [dict(seq.__dict__) for seq in ports(api).sequence.load_all(tags=list(tags))]
    except ApiEmptyResponse:
        return []


def analysis_records(api, tags: Sequence[str]) -> List[Dict[str, object]]:
    _, ApiEmptyResponse, _, _ = load_ibp_modules()
    try:
        return [dict(analyse.__dict__) for analyse in ports(api).zdna.load_all(tags=list(tags))]
    except ApiEmptyResponse:
        return []


def wait_for_jobs(api, jobs: Sequence[Dict[str, object]], job_type, poll_seconds: float) -> None:
    pending = {str(job["id"]): job for job in jobs}
    while pending:
        finished: List[str] = []
        for job_id, job in list(pending.items()):
            batch = ports(api).batch.get_batch_status(id=job_id, type=job_type)
            if batch.is_finished():
                finished.append(job_id)
            elif batch.is_failed():
                raise RuntimeError(f"IBP job failed: {job_id} {job.get('name', '')} {batch.exception or ''}")
        for job_id in finished:
            pending.pop(job_id, None)
        if pending:
            print(f"Waiting for {len(pending)} IBP job(s) ...")
            time.sleep(poll_seconds)


def write_record_fasta(record: FastaRecord, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f">{record.output_chrom}\n")
        for index in range(0, len(record.sequence), 80):
            handle.write(record.sequence[index : index + 80] + "\n")


def upload_sequences(api, records: Sequence[FastaRecord], args: argparse.Namespace) -> Dict[str, Dict[str, object]]:
    _, _, Types, normalize_name = load_ibp_modules()
    run_name = slug(args.run_name or args.fasta.stem)
    existing = {str(row.get("name")): row for row in sequence_records(api, ["zdna_fuzzy_detector"])} if args.reuse else {}
    jobs: List[Dict[str, object]] = []
    result: Dict[str, Dict[str, object]] = {}
    upload_dir = args.work_dir / "uploads" / run_name
    for record in records:
        seq_name = normalize_name(name=f"{run_name}_{slug(record.output_chrom)}_{record.coordinate_offset}_{len(record.sequence)}")
        if args.reuse and seq_name in existing:
            result[record.output_chrom] = existing[seq_name]
            print(f"Reusing sequence {seq_name}")
            continue
        fasta_path = upload_dir / f"{seq_name}.fa"
        write_record_fasta(record, fasta_path)
        seq = ports(api).sequence.create_file_sequence(
            circular=False,
            path=str(fasta_path.resolve()),
            name=seq_name,
            tags=["zdna_fuzzy_detector", run_name, record.output_chrom],
            nucleic_type="DNA",
            format="FASTA",
        )
        row = dict(seq.__dict__)
        result[record.output_chrom] = row
        jobs.append({"id": seq.id, "name": seq_name})
        print(f"Submitted sequence upload {seq_name}: {seq.id}")
    if jobs:
        wait_for_jobs(api, jobs, Types.SEQUENCE, args.poll_seconds)
    return result


def config_strictness(config: Dict[str, object]) -> Number:
    threshold = min(1.0, parse_float(config["threshold"]) / 70.0)
    size = min(1.0, parse_float(config["min_sequence_size"]) / 16.0)
    model_bonus = 1.0 if str(config["model"]) == "model2" else 0.85
    return 0.45 * threshold + 0.35 * size + 0.20 * model_bonus


def analyse_matches(row: Dict[str, object], config: Dict[str, object], sequence_id: str, tags: Sequence[str]) -> bool:
    row_tags = set(object_tags(row))
    return (
        str(row.get("sequence_id")) == str(sequence_id)
        and set(tags).issubset(row_tags)
        and str(row.get("model") or config["model"]) == str(config["model"])
        and parse_int(row.get("min_sequence_size")) == int(config["min_sequence_size"])
        and math.isclose(parse_float(row.get("min_score_percentage")), float(config["threshold"]), abs_tol=1e-6)
    )


def run_hunter_grid(api, sequences: Dict[str, Dict[str, object]], args: argparse.Namespace) -> Dict[Tuple[str, str], str]:
    _, _, Types, _ = load_ibp_modules()
    run_name = slug(args.run_name or args.fasta.stem)
    analysis_tag = f"zdna_fuzzy_{run_name}"
    existing = analysis_records(api, [analysis_tag]) if args.reuse else []
    jobs: List[Dict[str, object]] = []
    ids: Dict[Tuple[str, str], str] = {}
    for chrom, sequence in sequences.items():
        for config in HUNTER_CONFIGS:
            tags = [analysis_tag, str(config["config_id"]), chrom]
            match = next((row for row in existing if analyse_matches(row, config, str(sequence["id"]), tags)), None)
            if match is not None:
                analysis_id = str(match.get("id") or match.get("analysis_id"))
                print(f"Reusing analysis {chrom}/{config['config_id']}: {analysis_id}")
            else:
                analyse = ports(api).zdna.create_analyse(
                    id=str(sequence["id"]),
                    tags=tags,
                    min_sequence_size=int(config["min_sequence_size"]),
                    model=str(config["model"]),
                    GC_score=float(config["score_gc"]),
                    GTAC_score=float(config["score_gtac"]),
                    AT_score=float(config["score_at"]),
                    oth_score=float(config["score_oth"]),
                    min_score_percentage=float(config["threshold"]),
                )
                analysis_id = str(analyse.id)
                jobs.append({"id": analysis_id, "name": f"{chrom}/{config['config_id']}"})
                print(f"Submitted analysis {chrom}/{config['config_id']}: {analysis_id}")
            ids[(chrom, str(config["config_id"]))] = analysis_id
    if jobs:
        wait_for_jobs(api, jobs, Types.ZDNA, args.poll_seconds)
    return ids


def export_analysis(api, analysis_id: str) -> str:
    _, ApiEmptyResponse, _, _ = load_ibp_modules()
    try:
        return ports(api).zdna.export_csv(id=analysis_id) or ""
    except ApiEmptyResponse:
        return ""


def result_interval(row: Dict[str, str], coordinate_base: int, offset: int, chrom: str, source: str) -> Optional[Interval]:
    """Convert an API CSV row with inclusive end coordinates to half-open form."""
    fields = list(row.keys())
    start_col = detect_column(fields, ["start", "position", "pos", "from", "begin", "base_start", "start_position"])
    end_col = detect_column(fields, ["end", "stop", "to", "base_end", "end_position"])
    length_col = detect_column(fields, ["length", "size", "sequence_size", "window_length"])
    score_col = detect_column(fields, ["scoreperc", "score_percentage", "score", "value", "zdna_score", "z_score", "max_score"])
    if start_col is None:
        return None
    start = parse_int(row.get(start_col)) - coordinate_base + offset
    if end_col:
        end = parse_int(row.get(end_col)) - coordinate_base + offset
    elif length_col:
        end = start + max(1, parse_int(row.get(length_col))) - 1
    else:
        end = start
    if end < start:
        start, end = end, start
    score = parse_float(row.get(score_col), 0.0) if score_col else 0.0
    return Interval(chrom=chrom, start=start, end=end + 1, score=score, source=source)


def load_exported_intervals(
    api,
    analysis_ids: Dict[Tuple[str, str], str],
    records: Sequence[FastaRecord],
    args: argparse.Namespace,
) -> Dict[str, List[Interval]]:
    raw_dir = args.work_dir / "raw_exports" / slug(args.run_name or args.fasta.stem)
    raw_dir.mkdir(parents=True, exist_ok=True)
    by_record = {record.output_chrom: record for record in records}
    intervals_by_chrom: Dict[str, List[Interval]] = {}
    coordinate_base = int(args.api_coordinate_base)
    for (chrom, config_id), analysis_id in analysis_ids.items():
        raw_path = raw_dir / slug(chrom) / f"{config_id}.csv"
        raw_path.parent.mkdir(parents=True, exist_ok=True)
        if args.reuse and raw_path.exists():
            text = raw_path.read_text(encoding="utf-8", errors="replace")
            print(f"Reusing raw export {raw_path}")
        else:
            text = export_analysis(api, analysis_id)
            raw_path.write_text(text, encoding="utf-8")
        record = by_record[chrom]
        parsed = []
        for row in read_delimited_text(text):
            interval = result_interval(row, coordinate_base, record.coordinate_offset, chrom, config_id)
            if interval is not None:
                parsed.append(interval)
        intervals_by_chrom.setdefault(chrom, []).extend(parsed)
        print(f"Exported {len(parsed)} interval(s) from {chrom}/{config_id}")
    for intervals in intervals_by_chrom.values():
        intervals.sort(key=lambda item: (item.start, item.end))
    return intervals_by_chrom


def interval_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def parse_track_spec(spec: str, kind: str) -> Tuple[str, Path]:
    if "=" in spec:
        name, path_text = spec.split("=", 1)
        path = Path(path_text)
    else:
        path = Path(spec)
        name = path.stem
    name = slug(name)
    if not name:
        raise SystemExit(f"Invalid {kind} track name in: {spec}")
    if not path.exists():
        raise SystemExit(f"{kind.capitalize()} BED track not found: {path}")
    return name, path


def load_bed_intervals(path: Path) -> Dict[str, List[Tuple[int, int]]]:
    intervals: Dict[str, List[Tuple[int, int]]] = {}
    with path.open(encoding="utf-8-sig", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track ") or line.startswith("browser "):
                continue
            parts = re.split(r"[\t ]+", line)
            if len(parts) < 3:
                continue
            try:
                start = int(float(parts[1]))
                end = int(float(parts[2]))
            except ValueError:
                continue
            if end <= start:
                continue
            intervals.setdefault(canonical_chrom(parts[0]), []).append((start, end))
    for values in intervals.values():
        values.sort()
    return intervals


def load_annotation_tracks(specs: Sequence[str], kind: str) -> List[AnnotationTrack]:
    tracks: List[AnnotationTrack] = []
    seen: set[str] = set()
    for spec in specs:
        name, path = parse_track_spec(spec, kind)
        if name in seen:
            raise SystemExit(f"Duplicate {kind} track name: {name}")
        seen.add(name)
        tracks.append(AnnotationTrack(kind=kind, name=name, intervals_by_chrom=load_bed_intervals(path)))
    return tracks


def overlap_with_track(track: AnnotationTrack, chrom: str, start: int, end: int) -> Tuple[int, Number]:
    intervals = track.intervals_by_chrom.get(canonical_chrom(chrom), [])
    if not intervals:
        return 0, 0.0
    overlap_bp = 0
    for interval_start, interval_end in intervals:
        if interval_end <= start:
            continue
        if interval_start >= end:
            break
        overlap_bp += interval_overlap(start, end, interval_start, interval_end)
    return overlap_bp, clamp(100.0 * overlap_bp / max(1, end - start))


def tss_proximity_bin(distance_bp: int, promoter_window: int) -> int:
    if distance_bp <= 250:
        return 3
    if distance_bp <= promoter_window:
        return 2
    if distance_bp <= 5000:
        return 1
    return 0


def tss_context_label(distance_bp: int, promoter_window: int) -> str:
    if distance_bp == 0:
        return "overlaps_tss"
    if distance_bp <= 250:
        return "tss_250bp"
    if distance_bp <= promoter_window:
        return f"tss_{promoter_window}bp"
    if distance_bp <= 5000:
        return "tss_5kb"
    if distance_bp <= 20000:
        return "tss_20kb"
    return "distal"


def merge_candidate_components(intervals: Sequence[Interval]) -> List[Tuple[int, int, List[Interval]]]:
    """Merge overlapping half-open hits while retaining their source intervals."""
    if not intervals:
        return []
    sorted_intervals = sorted(intervals, key=lambda item: (item.start, item.end))
    components: List[Tuple[int, int, List[Interval]]] = []
    current_start = sorted_intervals[0].start
    current_end = sorted_intervals[0].end
    current_intervals = [sorted_intervals[0]]
    for interval in sorted_intervals[1:]:
        if interval.start <= current_end:
            current_end = max(current_end, interval.end)
            current_intervals.append(interval)
        else:
            components.append((current_start, current_end, current_intervals))
            current_start = interval.start
            current_end = interval.end
            current_intervals = [interval]
    components.append((current_start, current_end, current_intervals))
    return components


def component_stats_by_source(
    intervals: Sequence[Interval],
    start: int,
    end: int,
    sources: Sequence[str],
) -> Dict[str, Tuple[int, Number, int]]:
    stats: Dict[str, List[Number]] = {source: [0.0, 0.0, 0.0] for source in sources}
    for interval in intervals:
        if interval.source not in stats:
            continue
        overlap = interval_overlap(start, end, interval.start, interval.end)
        if overlap <= 0:
            continue
        current = stats[interval.source]
        current[0] += overlap
        current[1] = max(current[1], interval.score)
        current[2] = max(current[2], float(interval.length))
    return {source: (int(values[0]), values[1], int(values[2])) for source, values in stats.items()}


def piecewise_tss_score(distance_bp: Number) -> Number:
    distance = abs(distance_bp)
    if distance <= 250:
        return 100.0
    if distance <= 1000:
        return 90.0 - (distance - 250.0) * (20.0 / 750.0)
    if distance <= 5000:
        return 70.0 - (distance - 1000.0) * (35.0 / 4000.0)
    if distance <= 20000:
        return 35.0 - (distance - 5000.0) * (25.0 / 15000.0)
    return 0.0


def trapezoid_score(value: Number, a: Number, b: Number, c: Number, d: Number) -> Number:
    return 100.0 * FuzzyTerm("tmp", (a, b, c, d)).membership(value)


def degree(terms: Dict[str, FuzzyTerm], term: str, value: Number) -> Number:
    return terms[term].membership(value)


def components(features: Dict[str, Number], params: TunedParams) -> Dict[str, Number]:
    signal = clamp(
        params.signal_weighted * features["weighted"]
        + params.signal_zdna * features["zdna"]
        + params.signal_vote * features["hunter_vote"]
        + params.signal_model1 * features["model1"]
        + params.signal_model2 * features["model2"]
        + params.signal_all5_bonus * (features["all5"] / 100.0)
    )
    context = clamp(
        params.context_tss * features["tss"]
        + params.context_promoter * features["promoter"]
        + params.context_marks * features["marks"]
    )
    evidence = clamp(
        params.evidence_source * features["source_support"]
        + params.evidence_model1 * features["model1"]
        + params.evidence_model2 * features["model2"]
        + params.evidence_balance * features["balance"]
        + params.evidence_all5_bonus * (features["all5"] / 100.0)
    )
    bias = clamp(
        params.bias_repeat * features["repeat"]
        + params.bias_disagreement * features["disagreement"]
        + params.bias_model2_only * features["model2_only"]
        + params.bias_no_hit * features["no_hit"]
    )
    return {
        "signal": signal,
        "context": context,
        "evidence": evidence,
        "bias": bias,
        "feasibility": features["feasibility"],
        "balance": features["balance"],
    }


def fuzzy_rules(comp: Dict[str, Number]) -> Tuple[Number, List[Dict[str, object]]]:
    signal = comp["signal"]
    context = comp["context"]
    evidence = comp["evidence"]
    bias = comp["bias"]
    feasibility = comp["feasibility"]
    balance = comp["balance"]
    high_or_vh_evidence = max(degree(TERMS, "high", evidence), degree(TERMS, "very_high", evidence))
    low_or_vl_bias = max(degree(BIAS_TERMS, "very_low", bias), degree(BIAS_TERMS, "low", bias))
    med_or_high_context = max(degree(TERMS, "medium", context), degree(TERMS, "high", context), degree(TERMS, "very_high", context))
    low_or_med_evidence = max(degree(TERMS, "low", evidence), degree(TERMS, "medium", evidence))
    high_or_vh_bias = max(degree(BIAS_TERMS, "high", bias), degree(BIAS_TERMS, "very_high", bias))
    rules = [
        ("AF1", min(degree(TERMS, "very_high", signal), degree(TERMS, "very_high", evidence), degree(BIAS_TERMS, "very_low", bias)), 98.0),
        ("AF2", min(degree(TERMS, "very_high", signal), high_or_vh_evidence, low_or_vl_bias), 93.0),
        ("AF3", min(degree(TERMS, "high", signal), degree(TERMS, "very_high", evidence), med_or_high_context, low_or_vl_bias), 89.0),
        ("AF4", min(degree(TERMS, "high", signal), degree(TERMS, "high", evidence), max(degree(BIAS_TERMS, "low", bias), degree(BIAS_TERMS, "medium", bias))), 82.0),
        ("AF5", min(degree(TERMS, "medium", signal), high_or_vh_evidence, degree(TERMS, "high", context)), 72.0),
        ("AF6", min(degree(TERMS, "high", signal), degree(TERMS, "medium", evidence), degree(TERMS, "low", context)), 64.0),
        ("AF7", min(degree(TERMS, "medium", signal), degree(TERMS, "medium", evidence), degree(BIAS_TERMS, "low", bias)), 58.0),
        ("AF8", min(degree(TERMS, "very_low", signal), degree(TERMS, "low", evidence)), 12.0),
        ("AF9", min(high_or_vh_bias, low_or_med_evidence), 22.0),
        ("AF10", min(degree(TERMS, "high", signal), degree(TERMS, "low", balance)), 48.0),
        ("AF11", min(degree(TERMS, "very_high", signal), degree(TERMS, "very_high", balance), degree(TERMS, "medium", feasibility)), 91.0),
        ("AF12", min(degree(TERMS, "low", signal), high_or_vh_bias), 14.0),
    ]
    numerator = 0.0
    denominator = 0.0
    active: List[Dict[str, object]] = []
    for name, activation, consequent in rules:
        if activation <= 0.0:
            continue
        numerator += activation * consequent
        denominator += activation
        active.append({"rule": name, "activation": activation, "consequent": consequent})
    active.sort(key=lambda row: float(row["activation"]), reverse=True)
    return (50.0 if denominator == 0.0 else numerator / denominator), active


def model_specification() -> Dict[str, object]:
    """Return a machine-readable description of every fixed model element."""
    params = TunedParams()
    return {
        "name": "ZDNA-Fuzzy Hunter",
        "model_version": MODEL_VERSION,
        "coordinate_convention": "0-based, half-open [start, end)",
        "tuning_scope": (
            "Weights and the primary threshold were selected on the Shin human benchmark "
            "by deterministic random search guided by five stratified folds. The same folds "
            "contributed to hyperparameter selection, so the reported CV estimate is internal "
            "tuning-aware validation rather than nested or independent validation."
        ),
        "hunter_presets": HUNTER_PRESETS,
        "membership_functions": {
            "main_components": {name: [term.a, term.b, term.c, term.d] for name, term in TERMS.items()},
            "bias": {name: [term.a, term.b, term.c, term.d] for name, term in BIAS_TERMS.items()},
            "shape": "trapezoidal [a,b,c,d]",
        },
        "component_and_final_parameters": dict(params.__dict__),
        "rules": FUZZY_RULE_SPECS,
        "final_score": {
            "base": "offset + final_signal*signal + final_evidence*evidence + final_context*context + final_feasibility*feasibility - final_bias*bias",
            "mixture": "(1-rule_mix)*base + rule_mix*Takagi-Sugeno rule score",
            "rule_fraction": params.rule_mix,
            "weighted_component_fraction": 1.0 - params.rule_mix,
            "primary_threshold": FUZZY_THRESHOLD,
            "grid_hit_gate": True,
        },
        "operating_point_filters": {
            "moderate": {
                "exclude_if_tss_distance_lte_bp": MODERATE_TSS_BP,
                "and_score_lte": MODERATE_MAX_SCORE,
            },
            "strict": {
                "exclude_if_tss_distance_lte_bp": STRICT_TSS_BP,
                "and_signal_lte": STRICT_MAX_SIGNAL,
                "and_max_overlap_pct_lte": STRICT_MAX_OVERLAP_PCT,
            },
        },
    }


def emit_model_specification(path: Optional[Path]) -> None:
    """Print the model specification and optionally persist the same JSON."""
    text = json.dumps(model_specification(), indent=2)
    if path is not None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text + "\n", encoding="utf-8")
        print(f"Wrote model specification to {path}")
    else:
        print(text)


def score_features(features: Dict[str, Number], params: TunedParams) -> Tuple[Number, Dict[str, Number], List[Dict[str, object]]]:
    comp = components(features, params)
    rule_score, active_rules = fuzzy_rules(comp)
    base = clamp(
        params.offset
        + params.final_signal * comp["signal"]
        + params.final_evidence * comp["evidence"]
        + params.final_context * comp["context"]
        + params.final_feasibility * comp["feasibility"]
        - params.final_bias * comp["bias"]
    )
    score = clamp((1.0 - params.rule_mix) * base + params.rule_mix * rule_score)
    return score, comp, active_rules


def classify(mode: str, score: Number, comp: Dict[str, Number], tss_distance: int, max_overlap_pct: Number) -> bool:
    if score < FUZZY_THRESHOLD:
        return False
    if mode == "moderate" and tss_distance <= MODERATE_TSS_BP and score <= MODERATE_MAX_SCORE:
        return False
    if (
        mode == "strict"
        and tss_distance <= STRICT_TSS_BP
        and comp["signal"] <= STRICT_MAX_SIGNAL
        and max_overlap_pct <= STRICT_MAX_OVERLAP_PCT
    ):
        return False
    return True


def features_from_table_row(row: Dict[str, str]) -> Dict[str, Number]:
    has_grid = "grid_hit_fraction" in row
    if has_grid:
        weighted = 100.0 * parse_float(row.get("grid_weighted_hit_fraction"))
        hunter_vote = 100.0 * parse_float(row.get("grid_hit_fraction"))
        model1 = 100.0 * parse_float(row.get("grid_model1_hit_fraction"))
        model2 = 100.0 * parse_float(row.get("grid_model2_hit_fraction"))
        zdna = max(100.0 * parse_float(row.get("zdna_score")), parse_float(row.get("grid_max_raw_score")))
        strict = 100.0 * parse_float(row.get("grid_strict_hit_fraction"))
        consensus = 100.0 * parse_float(row.get("grid_model_consensus"))
        detection_class = row.get("grid_detection_class", "")
    else:
        selected = set(filter(None, row.get("selected_sources", "").split(";")))
        weighted = 100.0 * parse_float(row.get("weighted_hunter_signal"))
        hunter_vote = 100.0 * parse_float(row.get("hunter_vote_fraction"))
        model1 = 100.0 * parse_float(row.get("model1_vote_fraction"))
        model2 = 100.0 * parse_float(row.get("model2_vote_fraction"))
        zdna = 100.0 * parse_float(row.get("zdna_score"))
        strict = 100.0 if len(selected) >= 4 else 0.0
        consensus = 100.0 if len(selected) >= 5 else max(0.0, 100.0 - abs(model1 - model2))
        detection_class = ""
    tss_distance = parse_float(row.get("tss_distance_bp"), 50000.0)
    marks = clamp(
        100.0
        * parse_float(first_present(row, ["tss_proximity_bin", "cpx_context_count", "context_support_bin", "regulatory_marks"]))
        / 3.0
    )
    promoter = 100.0 if parse_bool(first_present(row, ["tss_proximal", "promoter_overlap"])) else 0.0
    repeat = clamp(4.0 * parse_float(first_present(row, ["heuristic_bias_score", "repeat_overlap_pct"])))
    disagreement = clamp(100.0 - consensus)
    all5 = strict
    model2_only = 100.0 if model2 > 0.0 and model1 == 0.0 else 0.0
    no_hit = 100.0 if weighted <= 0.0 else 0.0
    cross_tool = 100.0 if row.get("evidence_type") == "cross_tool" else 0.0
    predicted_motif = 100.0 if row.get("evidence_type") == "predicted_motif" else 0.0
    if has_grid:
        class_support = {
            "strict_consensus": 100.0,
            "stable_parameter_hit": 74.0,
            "permissive_only": 42.0,
            "no_grid_hit": 8.0,
        }.get(detection_class, 20.0)
        source_support = max(class_support, 0.70 * cross_tool + 0.30 * predicted_motif)
    else:
        source_support = 100.0 if cross_tool else 58.0 if predicted_motif else 14.0
    length = parse_float(row.get("length_bp"), 384.0)
    primer = parse_float(first_present(row, ["candidate_uniqueness_score", "primer_uniqueness"], 70.0), 70.0)
    feasibility = clamp(0.65 * trapezoid_score(length, 8.0, 8.0, 400.0, 450.0) + 0.35 * primer)
    return {
        "weighted": weighted,
        "hunter_vote": hunter_vote,
        "model1": model1,
        "model2": model2,
        "zdna": zdna,
        "tss": piecewise_tss_score(tss_distance),
        "marks": marks,
        "promoter": promoter,
        "repeat": repeat,
        "disagreement": disagreement,
        "balance": consensus,
        "all5": all5,
        "model2_only": model2_only,
        "no_hit": no_hit,
        "source_support": source_support,
        "feasibility": feasibility,
    }


def score_feature_table(args: argparse.Namespace) -> None:
    assert args.feature_table is not None
    assert args.output is not None
    started = time.perf_counter()
    rows = read_delimited_text(args.feature_table.read_text(encoding="utf-8-sig", errors="replace"))
    params = TunedParams()
    min_score = effective_min_score(args)
    scored_rows: List[Dict[str, object]] = []
    for row in rows:
        features = features_from_table_row(row)
        score, comp, active_rules = score_features(features, params)
        if args.require_grid_hit and parse_float(row.get("grid_hit_count")) <= 0.0:
            score = -1.0
        tss_distance = int(parse_float(row.get("tss_distance_bp"), 50000.0))
        max_overlap_pct = parse_float(row.get("grid_max_overlap_pct"))
        prediction = classify(args.mode, score, comp, tss_distance, max_overlap_pct) and score >= min_score
        out = rename_derived_feature_columns(row)
        out.update(
            {
                "preset": args.preset,
                "mode": args.mode,
                "min_score": f"{min_score:.4f}",
                "prediction": int(prediction),
                "fuzzy_score": f"{score:.4f}",
                "component_signal": f"{comp['signal']:.4f}",
                "component_context": f"{comp['context']:.4f}",
                "component_evidence": f"{comp['evidence']:.4f}",
                "component_bias": f"{comp['bias']:.4f}",
                "component_feasibility": f"{comp['feasibility']:.4f}",
                "component_balance": f"{comp['balance']:.4f}",
                "top_rules": ";".join(f"{item['rule']}:{float(item['activation']):.3f}" for item in active_rules[:3]),
            }
        )
        scored_rows.append(out)
    scoring_finished = time.perf_counter()
    written_rows = scored_rows if args.include_all else [row for row in scored_rows if int(row["prediction"]) == 1]
    out_format = infer_format(args.output, args.format)
    if out_format == "csv":
        columns = list(scored_rows[0].keys()) if scored_rows else []
        write_csv(written_rows, args.output, columns)
    else:
        bed_rows = []
        for row in written_rows:
            chrom = row.get("chrom") or row.get("chromosome") or row.get("context_chromosome")
            start = row.get("start") or row.get("shin_window_start") or row.get("context_window_start")
            end = row.get("end") or row.get("shin_window_end") or row.get("context_window_end")
            bed_rows.append({"chrom": chrom, "start": start, "end": end, "fuzzy_score": row["fuzzy_score"]})
        write_bedgraph(bed_rows, args.output)
    output_finished = time.perf_counter()
    timings = {
        "feature_table_load_and_scoring": round(scoring_finished - started, 6),
        "output_write": round(output_finished - scoring_finished, 6),
        "total_before_summary": round(output_finished - started, 6),
    }
    summary_data = summary(
        scored_rows,
        written_rows,
        args.mode,
        args.preset,
        min_score,
        args.hunter_config,
        backend="feature-table",
        timings=timings,
    )
    summary_path = args.summary_output or args.output.with_suffix(args.output.suffix + ".summary.json")
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary_data, indent=2), encoding="utf-8")
    print(f"Scored {len(scored_rows)} feature row(s).")
    print(f"Wrote {len(written_rows)} row(s) to {args.output}")
    print(f"Wrote summary to {summary_path}")


def candidate_rows(
    intervals_by_chrom: Dict[str, List[Interval]],
    tss_by_chrom: Dict[str, List[int]],
    annotation_tracks: Sequence[AnnotationTrack],
    args: argparse.Namespace,
) -> List[Dict[str, object]]:
    params = TunedParams()
    min_score = effective_min_score(args)
    total_weight = sum(config_strictness(config) for config in HUNTER_CONFIGS)
    strict_sources = {config["config_id"] for config in HUNTER_CONFIGS if config["threshold"] >= 60 or config["min_sequence_size"] >= 10}
    config_ids = [str(config["config_id"]) for config in HUNTER_CONFIGS]
    rows: List[Dict[str, object]] = []
    for chrom in sorted(intervals_by_chrom):
        intervals = intervals_by_chrom[chrom]
        components = merge_candidate_components(intervals)
        print(
            f"Building candidates for {chrom}: {len(intervals)} interval(s) -> {len(components)} merged candidate(s)",
            flush=True,
        )
        tss_positions = tss_by_chrom.get(canonical_chrom(chrom), [])
        for index, (start, end, component_intervals) in enumerate(components, start=1):
            if index % 100000 == 0:
                print(f"  scored {index} candidate(s) for {chrom} ...", flush=True)
            length = max(1, end - start)
            hit_count = 0
            strict_hit_count = 0
            weighted_hit = 0.0
            overlap_values: List[Number] = []
            score_values: List[Number] = []
            best_interval_length = 0
            flags: Dict[str, int] = {}
            overlap_by_config: Dict[str, Number] = {}
            source_stats = component_stats_by_source(component_intervals, start, end, config_ids)
            for config in HUNTER_CONFIGS:
                config_id = str(config["config_id"])
                overlap_bp, max_score, best_length = source_stats[config_id]
                hit = overlap_bp > 0
                flags[config_id] = int(hit)
                overlap_pct = 100.0 * overlap_bp / length
                overlap_by_config[config_id] = overlap_pct
                if hit:
                    hit_count += 1
                    weighted_hit += config_strictness(config)
                    overlap_values.append(overlap_pct)
                    score_values.append(max_score)
                    best_interval_length = max(best_interval_length, best_length)
                    if config_id in strict_sources:
                        strict_hit_count += 1
            hit_fraction = hit_count / max(1, len(HUNTER_CONFIGS))
            weighted_fraction = weighted_hit / max(total_weight, 1e-9)
            strict_fraction = strict_hit_count / max(1, len(HUNTER_CONFIGS))
            max_overlap = max(overlap_values) if overlap_values else 0.0
            mean_overlap = sum(overlap_values) / len(overlap_values) if overlap_values else 0.0
            max_score = max(score_values) if score_values else 0.0
            if hit_count == 0:
                detection_class = "no_grid_hit"
            elif strict_hit_count > 0 and hit_count == len(HUNTER_CONFIGS):
                detection_class = "strict_consensus"
            elif weighted_fraction >= 0.25:
                detection_class = "stable_parameter_hit"
            else:
                detection_class = "permissive_only"
            tss_distance = nearest_distance(tss_positions, start, end)
            tss_distance = 50000 if tss_distance is None else int(tss_distance)
            proximity_bin = tss_proximity_bin(tss_distance, args.promoter_window)
            tss_proximal = int(tss_distance <= args.promoter_window)
            features = {
                "weighted": 100.0 * weighted_fraction,
                "hunter_vote": 100.0 * hit_fraction,
                "model1": 0.0,
                "model2": 100.0 * hit_fraction,
                "zdna": max_score,
                "tss": piecewise_tss_score(tss_distance),
                "marks": clamp(100.0 * proximity_bin / 3.0),
                "promoter": 100.0 if tss_proximal else 0.0,
                "repeat": 0.0,
                "disagreement": 100.0,
                "balance": 0.0,
                "all5": 100.0 * strict_fraction,
                "model2_only": 100.0 if hit_count > 0 else 0.0,
                "no_hit": 100.0 if hit_count == 0 else 0.0,
                "source_support": {"strict_consensus": 100.0, "stable_parameter_hit": 74.0, "permissive_only": 42.0, "no_grid_hit": 8.0}[detection_class],
                "feasibility": clamp(0.65 * trapezoid_score(length, 8.0, 8.0, 400.0, 450.0) + 0.35 * 70.0),
            }
            score, comp, active_rules = score_features(features, params)
            prediction = classify(args.mode, score, comp, tss_distance, max_overlap) and score >= min_score
            row = {
                "chrom": chrom,
                "start": start,
                "end": end,
                "candidate_id": f"{slug(chrom)}_{index}",
                "preset": args.preset,
                "mode": args.mode,
                "min_score": round(min_score, 4),
                "prediction": int(prediction),
                "fuzzy_score": round(score, 4),
                "component_signal": round(comp["signal"], 4),
                "component_context": round(comp["context"], 4),
                "component_evidence": round(comp["evidence"], 4),
                "component_bias": round(comp["bias"], 4),
                "component_feasibility": round(comp["feasibility"], 4),
                "component_balance": round(comp["balance"], 4),
                "length_bp": length,
                "tss_distance_bp": tss_distance,
                "tss_proximal": tss_proximal,
                "tss_proximity_bin": proximity_bin,
                "tss_context_label": tss_context_label(tss_distance, args.promoter_window),
                "grid_hit_count": hit_count,
                "grid_hit_fraction": round(hit_fraction, 6),
                "grid_weighted_hit_fraction": round(weighted_fraction, 6),
                "grid_model1_hit_fraction": 0.0,
                "grid_model2_hit_fraction": round(hit_fraction, 6),
                "grid_model_consensus": 0.0,
                "grid_strict_hit_fraction": round(strict_fraction, 6),
                "grid_max_overlap_pct": round(max_overlap, 4),
                "grid_mean_overlap_pct": round(mean_overlap, 4),
                "grid_max_raw_score": round(max_score, 4),
                "grid_best_interval_length": best_interval_length,
                "grid_detection_class": detection_class,
                "top_rules": ";".join(f"{item['rule']}:{float(item['activation']):.3f}" for item in active_rules[:3]),
            }
            for config in HUNTER_CONFIGS:
                config_id = str(config["config_id"])
                row[f"grid_hit__{config_id}"] = flags[config_id]
                row[f"grid_overlap_pct__{config_id}"] = round(overlap_by_config[config_id], 4)
            for track in annotation_tracks:
                overlap_bp, overlap_pct = overlap_with_track(track, chrom, start, end)
                prefix = "epigenomic" if track.kind == "epigenomic" else "annotation"
                row[f"{prefix}_{track.name}_hit"] = int(overlap_bp > 0)
                row[f"{prefix}_{track.name}_overlap_pct"] = round(overlap_pct, 4)
            rows.append(row)
        print(f"Finished {chrom}: scored {len(components)} candidate(s)", flush=True)
    rows.sort(key=lambda row: (str(row["chrom"]), int(row["start"]), int(row["end"])))
    return rows


def output_columns() -> List[str]:
    base = [
        "chrom",
        "start",
        "end",
        "candidate_id",
        "preset",
        "mode",
        "min_score",
        "prediction",
        "fuzzy_score",
        "component_signal",
        "component_context",
        "component_evidence",
        "component_bias",
        "component_feasibility",
        "component_balance",
        "length_bp",
        "tss_distance_bp",
        "tss_proximal",
        "tss_proximity_bin",
        "tss_context_label",
        "grid_hit_count",
        "grid_hit_fraction",
        "grid_weighted_hit_fraction",
        "grid_model1_hit_fraction",
        "grid_model2_hit_fraction",
        "grid_model_consensus",
        "grid_strict_hit_fraction",
        "grid_max_overlap_pct",
        "grid_mean_overlap_pct",
        "grid_max_raw_score",
        "grid_best_interval_length",
        "grid_detection_class",
    ]
    for config in HUNTER_CONFIGS:
        base.append(f"grid_hit__{config['config_id']}")
        base.append(f"grid_overlap_pct__{config['config_id']}")
    base.append("top_rules")
    return base


def write_bedgraph(rows: Sequence[Dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(f"{row['chrom']}\t{row['start']}\t{row['end']}\t{row['fuzzy_score']}\n")


def infer_format(path: Path, explicit: Optional[str]) -> str:
    if explicit:
        return explicit
    suffix = path.suffix.lower()
    if suffix in {".bedgraph", ".bdg"}:
        return "bedgraph"
    return "csv"


def summary(
    rows: Sequence[Dict[str, object]],
    written_rows: Sequence[Dict[str, object]],
    mode: str,
    preset: str,
    min_score: Number,
    hunter_config_path: Optional[Path],
    backend: str,
    timings: Optional[Dict[str, Number]] = None,
) -> Dict[str, object]:
    scores = [float(row["fuzzy_score"]) for row in rows]
    by_class: Dict[str, int] = {}
    by_chrom: Dict[str, int] = {}
    for row in rows:
        detection_class = str(row.get("grid_detection_class", "unknown"))
        chrom = str(row.get("chrom") or row.get("chromosome") or row.get("context_chromosome") or "unknown")
        by_class[detection_class] = by_class.get(detection_class, 0) + 1
        by_chrom[chrom] = by_chrom.get(chrom, 0) + 1
    return {
        "tool": "ZDNA-Fuzzy Hunter",
        "model_version": MODEL_VERSION,
        "hunter_backend": backend,
        "preset": preset,
        "hunter_config_file": str(hunter_config_path) if hunter_config_path else "",
        "mode": mode,
        "min_score": round(min_score, 4),
        "total_candidates": len(rows),
        "written_rows": len(written_rows),
        "selected_candidates": sum(1 for row in rows if int(row["prediction"]) == 1),
        "mean_fuzzy_score": round(statistics.mean(scores), 4) if scores else 0.0,
        "median_fuzzy_score": round(statistics.median(scores), 4) if scores else 0.0,
        "max_fuzzy_score": round(max(scores), 4) if scores else 0.0,
        "detection_class_counts": by_class,
        "chromosome_candidate_counts": by_chrom,
        "hunter_configs": HUNTER_CONFIGS,
        "fuzzy_threshold": FUZZY_THRESHOLD,
        "runtime": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "timings_seconds": timings or {},
        },
    }


def maybe_cleanup(args: argparse.Namespace) -> None:
    if args.keep_intermediate:
        return
    # Keep raw exports by default because they are useful for reproducibility.
    # Uploaded FASTA slices can contain genome data and are removed unless requested.
    upload_dir = args.work_dir / "uploads" / slug(args.run_name or args.fasta.stem)
    if not upload_dir.exists():
        return
    for path in sorted(upload_dir.glob("*.fa")):
        path.unlink()


def main() -> None:
    started = time.perf_counter()
    args = parse_args()
    apply_preset(args)
    if args.describe_model:
        emit_model_specification(args.model_spec_output)
        return
    if args.model_spec_output is not None:
        emit_model_specification(args.model_spec_output)
    if args.output is None:
        raise SystemExit("--output is required unless --describe-model is used.")
    if args.feature_table is not None:
        score_feature_table(args)
        return
    if args.fasta is None or args.tss is None:
        raise SystemExit("--fasta and --tss are required unless --score-feature-table is used.")
    args.work_dir.mkdir(parents=True, exist_ok=True)
    records = selected_fasta_records(args)
    tss = load_tss(args.tss, args.tss_coordinate_base)
    genomic_tracks = load_annotation_tracks(args.annotation_bed, "genomic")
    epigenomic_tracks = load_annotation_tracks(args.epigenomic_bed, "epigenomic")
    annotation_tracks = [*genomic_tracks, *epigenomic_tracks]
    print(f"Loaded {len(records)} FASTA record(s) and TSS annotations for {len(tss)} chromosome(s).")
    if genomic_tracks:
        print(f"Loaded {len(genomic_tracks)} genomic-context BED annotation track(s).")
    if epigenomic_tracks:
        print(f"Loaded {len(epigenomic_tracks)} epigenomic BED annotation track(s).")
    config_source = f"custom hunter config {args.hunter_config}" if args.hunter_config else f"preset {args.preset}"
    print(
        f"Using {config_source} with {len(HUNTER_CONFIGS)} Z-DNA Hunter config(s), "
        f"backend={args.hunter_backend}, mode={args.mode}, min_score={effective_min_score(args):.4f}."
    )
    input_finished = time.perf_counter()
    if args.hunter_backend == "local":
        intervals_by_chrom = run_local_hunter(records)
    else:
        api = make_api(args)
        sequences = upload_sequences(api, records, args)
        analysis_ids = run_hunter_grid(api, sequences, args)
        intervals_by_chrom = load_exported_intervals(api, analysis_ids, records, args)
    hunter_finished = time.perf_counter()
    rows = candidate_rows(intervals_by_chrom, tss, annotation_tracks, args)
    scoring_finished = time.perf_counter()
    written_rows = rows if args.include_all else [row for row in rows if int(row["prediction"]) == 1]
    out_format = infer_format(args.output, args.format)
    if out_format == "csv":
        write_csv(written_rows, args.output)
    else:
        write_bedgraph(written_rows, args.output)
    output_finished = time.perf_counter()
    timings = {
        "input_loading": round(input_finished - started, 6),
        "hunter": round(hunter_finished - input_finished, 6),
        "candidate_scoring": round(scoring_finished - hunter_finished, 6),
        "output_write": round(output_finished - scoring_finished, 6),
        "total_before_summary": round(output_finished - started, 6),
    }
    summary_data = summary(
        rows,
        written_rows,
        args.mode,
        args.preset,
        effective_min_score(args),
        args.hunter_config,
        backend=args.hunter_backend,
        timings=timings,
    )
    summary_path = args.summary_output or args.output.with_suffix(args.output.suffix + ".summary.json")
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary_data, indent=2), encoding="utf-8")
    if args.hunter_backend == "api":
        maybe_cleanup(args)
    print(f"Wrote {len(written_rows)} row(s) to {args.output}")
    print(f"Wrote summary to {summary_path}")


if __name__ == "__main__":
    main()
