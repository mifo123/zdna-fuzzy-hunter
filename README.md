# ZDNA-Fuzzy Hunter

ZDNA-Fuzzy Hunter is a command-line tool for fast, interpretable prioritization
of Z-DNA-forming regions. It combines two complementary Z-DNA Hunter scans with
a fixed fuzzy expert layer using sequence evidence and nearest-transcription
start site (TSS) context.

The default workflow is fully local and uses only the Python standard library.
An optional IBP API backend is retained for hosted Z-DNA Hunter analyses.

## What the tool does

The `genome-balanced` preset runs two Z-DNA Hunter model 2 configurations:

| ID | Minimum length | Threshold | Role |
|---|---:|---:|---|
| `m2_l8_t30` | 8 bp | 30% | permissive candidate recovery |
| `m2_l10_t60` | 10 bp | 60% | stringent confirmation |

The `shin-publication` preset changes the permissive arm to `m2_l6_t30` to
reproduce the Shin benchmark configuration. The `genome-strict` preset runs only
`m2_l10_t60`.

Overlapping scan hits are merged. For every candidate, the tool computes:

- per-configuration hit and overlap statistics;
- maximum raw Z-DNA Hunter score and candidate length;
- distance to the nearest TSS;
- sequence-signal, evidence, context, bias and feasibility components;
- the fuzzy score and the strongest activated rules;
- optional overlaps with user-supplied genomic or epigenomic BED tracks.

All exported genomic intervals use **0-based, half-open `[start, end)`
coordinates**. A region entered with `--region`, such as
`chr2L:1000-5000`, is deliberately interpreted as a **1-based closed** user
interval and converted internally.

## Local and API backends

The default local backend is a Python port of the Z-DNA Hunter state machine
used by the MENDELU/IBP backend. It runs the complete FASTA-to-candidate pipeline
without uploading sequences:

```text
FASTA + TSS -> local paired Z-DNA Hunter scans -> fuzzy layer -> CSV/bedGraph
```

The local implementation was regression-checked on the eight main Drosophila
dm6 chromosomes. It reproduced 482,286 merged candidates and 33,125 strict calls
from the stored IBP run. Coordinates and calls matched; 16 of 33,125 exported
scores differed only by 0.0001 because the API CSV stores rounded raw scores.

Use `--hunter-backend api` when a hosted IBP run is preferred. API use requires
credentials and the optional client dependency described below. Raw API exports
are retained for provenance.

## Fuzzy model and parameters

Five trapezoidal linguistic values are used for the main components: very low,
low, medium, high and very high. Bias uses a separate set of five trapezoids.
Twelve Takagi-Sugeno-style rules produce a rule score. The final score is:

```text
base = offset
       + 0.439786 * signal
       + 0.448789 * evidence
       + 0.010134 * context
       + 0.101291 * feasibility
       - 0.471126 * bias

final = 0.850764 * base + 0.149236 * rule_score
```

Thus, the rule output contributes **14.9236%**, and the weighted component score
contributes **85.0764%**, to the final mixture. The publication threshold is
19.157, with an additional gate requiring at least one upstream Hunter hit.

All 12 fuzzy rules are evaluated for every candidate in every operating mode.
The balanced, moderate and strict modes do not enable or disable rules; they
apply different decision filters only after the fuzzy score has been computed.
In rule descriptions, a slash or an `_or_` expression denotes logical OR. On
the released Shin feature table, AF5-AF7 and AF10-AF11 have zero activation
because no locus satisfies their antecedent memberships, not because a mode
switches them off.

Every component coefficient, membership breakpoint, operating-point filter and
rule is available in machine-readable JSON:

```bash
python zdna_fuzzy_detector.py --describe-model
python zdna_fuzzy_detector.py \
  --describe-model \
  --model-spec-output model_specification.json
```

The model specification is generated from the same constants used by the
scoring code, preventing documentation from drifting away from the executable
model.

## Parameter selection and validation scope

The fixed fuzzy weights were selected by deterministic random search on the
Shin human benchmark. Candidate weight sets were ranked using predictions from
five stratified folds, and each fold selected its decision threshold from its
training portion. However, the same five folds contributed to weight selection
and to the reported cross-validated estimate. The result is therefore an
**internal, tuning-aware validation**, not fully nested cross-validation and not
an independent estimate of generalization.

The reported Shin operating points are:

| Model / mode | Recall | Specificity | Balanced accuracy |
|---|---:|---:|---:|
| Z-DNABERT hg18 threshold 0.25, published external row | 0.877 | 0.890 | 0.880 |
| ZDNA-Fuzzy Hunter, balanced | 0.862 | 0.904 | 0.883 |
| ZDNA-Fuzzy Hunter, moderate-specificity | 0.843 | 0.945 | 0.894 |
| ZDNA-Fuzzy Hunter, strict-specificity | 0.833 | 0.973 | 0.903 |

Because Z-DNABERT was not rerun in the same pipeline, the near-equal balanced
accuracy values should be interpreted as comparable published operating points,
not evidence that one method outperforms the other.

### Audited Shin analysis reproduction

The public HG Shin FASTA cohort linked by the Z-DNABERT repository contains 391
records (316 positive and 75 negative). The reported analysis used 385 records
on the 24 primary hg18 chromosomes; six `*_random`-contig records were excluded.
The exact identifiers and reason are machine-readable in
[`shin_expected_metrics.json`](validation/data/shin_expected_metrics.json).

The repository includes those source FASTA files, the author-derived context
variables used during model fitting, and a standard-library-only reproduction
script. From a fresh clone, run:

```bash
python validation/reproduce_shin_benchmark.py
```

The command audits the 391-to-385 cohort accounting, converts the source FASTA
headers from 1-based closed to 0-based half-open coordinates, rebuilds both
Z-DNA Hunter feature sets locally, and then evaluates the balanced, moderate
and strict modes. It exits with an error if the rebuilt feature table or any
metric differs from the released expected values.

`shin_context_features.csv` is explicitly an author-generated intermediate,
not original Shin data and not independent validation. Its role and derivation
limits are documented in
[`validation/data/README.md`](validation/data/README.md).

### Independent U2OS check

The fixed Shin model was also checked against the reviewer-suggested human U2OS
Z-DNA ChIP-seq series **GSE290662**. The analysis used hg19 chromosome 22,
GENCODE v19 TSS annotation, two WT Z22-IP replicates (`GSM8818990` and
`GSM8818991`) and two matched input replicates (`GSM8818982` and `GSM8818983`).
No U2OS measurement was used to change a weight, rule or threshold.

Among 20,000 candidates sampled deterministically across fuzzy-score deciles,
10.05% of the highest-decile windows and 7.40% of the lowest-decile windows met
the replicate-consistent high-Z22-IP criterion (risk ratio 1.36; bootstrap 95%
CI for the rate difference 0.009-0.044). The continuous fuzzy-score versus
Z22-IP/input rank correlation was 0.0004, and the high-decile enrichment was not
observed for candidates more than 20 kb from a TSS. This is therefore modest,
TSS-associated independent support, not evidence of broad assay or species
transfer. Exact inputs, deterministic sampling and output fields are documented
in [`validation/README.md`](validation/README.md).

### Independent original-data rice check

We also used the original data associated with reviewer-specified
**doi:10.1111/pbi.14585**: rice ZIP-seq/CUT&Tag series **GSE252954**. The unchanged model generated
558,171 candidates on MSU/UGA Release 7 chromosome 1. In a deterministic
20,000-candidate sample, 16.8% of the highest-decile and 2.8% of the
lowest-decile windows met the prespecified replicate-consistent ZIP-seq
criterion (risk ratio 6.0; bootstrap 95% CI for the rate difference
0.121-0.158), although the complete decile pattern was non-monotonic.
Continuous fuzzy score correlated modestly with ZIP-seq
enrichment (Spearman 0.162) and weakly with CUT&Tag signal (0.065). Every
highest-decile sampled candidate was within 20 kb of a TSS, so this supports a
TSS-proximal cross-species signal but does not isolate sequence transfer from
the model's explicit TSS feature.

## Installation

Python **3.10 or newer** is required. The code has been tested with Python 3.10
and 3.12.

The local backend has no third-party runtime dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
```

For the optional IBP API backend:

```bash
pip install -r requirements-api.txt
export IBP_EMAIL="your.email@example.org"
export IBP_PASSWORD="your-password"
```

Credentials are read from `--email` / `--password` or from `IBP_EMAIL` /
`IBP_PASSWORD`; they are not written to output files.

## Input formats

### FASTA

Single-record and multi-record FASTA files are accepted:

```text
>chr2L
ACGT...
>chr2R
ACGT...
```

Use `--chromosomes chr2L,chr2R` to restrict records or
`--region chr2L:100000-200000` to scan a 1-based closed interval.

### TSS annotation

SGA, BED and headered CSV/TSV files are supported. SGA input is recognized
directly, for example:

```text
chr1  TSS  850984  +  1  SAMD11
```

BED input uses ordinary 0-based half-open coordinates:

```text
chr1  850983  850984  SAMD11  0  +
```

With the default `--tss-coordinate-base auto`, SGA point coordinates are
treated as 1-based and BED intervals as 0-based. For a negative-strand BED
feature, the TSS is the last covered base (`end - 1`). For ambiguous generic
point files, set `--tss-coordinate-base 0` or `1` explicitly.

### Optional annotation tracks

BED, narrowPeak and broadPeak tracks may be attached for output annotation.
They do not change the published fuzzy score unless the model is explicitly
retuned.

```bash
python zdna_fuzzy_detector.py \
  --fasta data/genome.fa \
  --tss data/tss.sga \
  --annotation-bed promoters=tracks/promoters.bed \
  --epigenomic-bed ATAC=tracks/atac_peaks.narrowPeak \
  --output results/zdna_annotated.csv
```

## Common commands

Fully local whole-genome analysis:

```bash
python zdna_fuzzy_detector.py \
  --fasta data/genome.fa \
  --tss data/tss.sga \
  --preset genome-balanced \
  --mode balanced \
  --output results/zdna_balanced.csv
```

Strict calls on selected chromosomes:

```bash
python zdna_fuzzy_detector.py \
  --fasta data/genome.fa \
  --tss data/tss.sga \
  --chromosomes chr2L,chr2R \
  --preset genome-balanced \
  --mode strict \
  --min-score 50 \
  --output results/dm6_strict.csv
```

Hosted IBP backend:

```bash
python zdna_fuzzy_detector.py \
  --hunter-backend api \
  --fasta data/genome.fa \
  --tss data/tss.sga \
  --preset genome-balanced \
  --mode balanced \
  --output results/zdna_api.csv
```

Custom Hunter grid:

```csv
config_id,model,min_sequence_size,threshold,score_gc,score_gtac,score_at,score_oth
custom_l8_t40,model2,8,40,2,1,0.5,0
custom_l10_t65,model2,10,65,2,1,0.5,0
```

```bash
python zdna_fuzzy_detector.py \
  --fasta data/genome.fa \
  --tss data/tss.sga \
  --hunter-config custom_hunter.csv \
  --mode balanced \
  --output results/custom.csv
```

Use unique `config_id` values. For API reruns with changed parameters, also use
a new `--run-name` or `--no-reuse`.

## Outputs

CSV output includes:

- `chrom`, `start`, `end`, `candidate_id`;
- `preset`, `mode`, `min_score`, `prediction`, `fuzzy_score`;
- component scores and nearest-TSS fields;
- grid hit, overlap, score and detection-class fields;
- optional `annotation_*` and `epigenomic_*` columns;
- top active fuzzy rules.

bedGraph output contains four columns:

```text
chrom  start  end  fuzzy_score
```

By default, only accepted candidates are written. `--include-all` exports every
merged candidate. Each run also writes a JSON summary containing the complete
Hunter configuration, model version, backend, counts, hardware/Python metadata
and stage-level timings.

### Feature-table reproducibility mode

Prepared candidate features can be rescored without running either Hunter
backend:

```bash
python zdna_fuzzy_detector.py \
  --score-feature-table validation/data/shin_analysis_features.csv \
  --preset shin-publication \
  --mode balanced \
  --include-all \
  --output results/shin_rescored.csv
```

`supplementary/Supplementary_Table_S1_Shin_per_locus.csv` is a compact
publication output, not a scoring input. The canonical reproduction command is
`python validation/reproduce_shin_benchmark.py`, which first rebuilds the
analysis table from the released source and context inputs.

## Runtime benchmark

One complete local run on the eight main dm6 chromosomes used a MacBook Air
with an Apple M2 (8 CPU cores), 16 GB RAM, macOS 26.6 and Python 3.12.13.
The reported benchmark used the `genome-balanced` preset
(`m2_l8_t30` plus `m2_l10_t60`), `strict` mode and `--min-score 50`:

```bash
python zdna_fuzzy_detector.py \
  --fasta data/dm6/dm6.fa \
  --tss data/dm6/dm6_refGene_tss.sga \
  --chromosomes chr2L,chr2R,chr3L,chr3R,chr4,chrM,chrX,chrY \
  --preset genome-balanced \
  --mode strict \
  --min-score 50 \
  --output results/dm6_main_chromosomes_strict.csv
```

| Stage | Time (s) |
|---|---:|
| FASTA/TSS loading | 1.875 |
| two local Z-DNA Hunter scans | 77.457 |
| merging and fuzzy scoring of 482,286 candidates | 8.975 |
| writing 33,125 strict calls | 0.470 |
| measured total before JSON summary | 88.776 |
| external wall-clock total | 89.580 |

Separately, rescoring and writing the existing 33,125-row dm6 feature table took
a median 0.91 s over five runs (range 0.90-0.96 s). Runtime depends on genome
size, number of Hunter configurations, storage and optional API/network latency;
therefore the tool does not claim that every whole-genome run completes within
seconds.

## Testing

Run the standard-library test suite with:

```bash
python -m unittest discover -s tests -v
```

Tests cover SGA/BED coordinate normalization, negative-strand BED TSS handling,
API inclusive-end conversion, local Hunter coordinates and model-specification
completeness. They also rerun all three Shin operating points from the released
feature table and verify the exact confusion matrices and metrics.

## Limitations

- Fuzzy weights and the primary threshold were developed on the human Shin
  benchmark. The GSE290662 and original-data GSE252954 checks are limited to one
  chromosome each, deterministic samples and fixed-window signal criteria; they
  are not second labelled classification benchmarks. The rice result is
  cross-species but fully TSS-proximal in its highest decile, so broad species
  and threshold transfer still require prospective testing.
- TSS proximity is associated with the Shin benchmark and can lift borderline
  candidates. Inspect component scores and use sequence-only evidence when
  comparing datasets with very different TSS distributions.
- Context and epigenomic BED overlaps are annotations unless a model is retuned
  to use them.
- Z-DNA Hunter, Z-Hunt II/Z-GENIE, ZSeeker, DeepZ and Z-DNABERT implement
  different scoring assumptions; cross-tool numerical scores are not directly
  interchangeable.

## License

ZDNA-Fuzzy Hunter is released under the [MIT License](LICENSE).

## Citation and preservation

Machine-readable citation metadata are provided in [`CITATION.cff`](CITATION.cff).
Version 1.1.0, including the exact code, validation inputs and supplementary
files, is preserved on Zenodo. Please cite the version-specific DOI:
[10.5281/zenodo.22279840](https://doi.org/10.5281/zenodo.22279840).
