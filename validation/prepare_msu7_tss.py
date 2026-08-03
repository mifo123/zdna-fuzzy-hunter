#!/usr/bin/env python3
"""Convert MSU/UGA Rice Genome Annotation Release 7 mRNAs to TSS CSV."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gff3", type=Path, required=True)
    parser.add_argument("--chromosome", default="Chr1")
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def parse_attributes(raw: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for item in raw.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attributes[key] = value
    return attributes


def main() -> None:
    args = parse_args()
    rows: list[dict[str, object]] = []
    seen: set[tuple[str, int, str, str]] = set()
    with gzip.open(args.gff3, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[0] != args.chromosome or fields[2] != "mRNA":
                continue
            start_1, end_1 = int(fields[3]), int(fields[4])
            strand = fields[6]
            tss_0 = start_1 - 1 if strand == "+" else end_1 - 1
            attrs = parse_attributes(fields[8])
            transcript_id = attrs.get("ID", "")
            gene_id = attrs.get("Parent", "")
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
                    "gene_id": gene_id,
                    "gene_name": gene_id,
                    "transcript_id": transcript_id,
                    "annotation": "MSU/UGA Rice Genome Annotation Release 7",
                }
            )
    if not rows:
        raise SystemExit(f"No mRNA records found for {args.chromosome}")
    rows.sort(key=lambda row: (int(row["start"]), str(row["strand"]), str(row["transcript_id"])))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows):,} TSS records to {args.output}")


if __name__ == "__main__":
    main()
