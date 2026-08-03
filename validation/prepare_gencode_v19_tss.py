#!/usr/bin/env python3
"""Convert transcript records in a gzipped GENCODE GTF to a TSS CSV."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--chromosome", default="chr22")
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def parse_attributes(raw: str) -> dict[str, str]:
    return {key: value for key, value in re.findall(r'(\S+)\s+"([^"]*)"', raw)}


def main() -> None:
    args = parse_args()
    rows: list[dict[str, object]] = []
    seen: set[tuple[str, int, str, str]] = set()
    with gzip.open(args.gtf, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[0] != args.chromosome or fields[2] != "transcript":
                continue
            start_1, end_1 = int(fields[3]), int(fields[4])
            strand = fields[6]
            tss_0 = start_1 - 1 if strand == "+" else end_1 - 1
            attrs = parse_attributes(fields[8])
            transcript_id = attrs.get("transcript_id", "")
            key = (fields[0], tss_0, strand, transcript_id)
            if key in seen:
                continue
            seen.add(key)
            rows.append(
                {
                    "chromosome": fields[0],
                    "start": tss_0,
                    "end": tss_0 + 1,
                    "strand": strand,
                    "gene_id": attrs.get("gene_id", ""),
                    "gene_name": attrs.get("gene_name", ""),
                    "transcript_id": transcript_id,
                    "annotation": "GENCODE v19",
                }
            )
    if not rows:
        raise SystemExit(f"No transcript records found for {args.chromosome}")
    rows.sort(key=lambda row: (int(row["start"]), str(row["strand"]), str(row["transcript_id"])))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows):,} TSS records to {args.output}")


if __name__ == "__main__":
    main()
