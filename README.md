# Z-DNA Fuzzy Detector

Command line tool for fast genome-scale prioritization of Z-DNA forming
regions. The pipeline runs selected Z-DNA Hunter analyses through the IBP API
and applies a tuned fuzzy expert layer using sequence signal and nearest-TSS
proximity context.

The tool is intended as a lightweight, interpretable layer above Z-DNA Hunter.
It accepts one chromosome, a multi-chromosome genome FASTA, or a selected
subregion and exports prioritized candidates as CSV or bedGraph.

## Method

The tool provides several Z-DNA Hunter presets. The default is
`genome-balanced`, intended for practical genome-scale use:

| ID | Model | Minimal length | Threshold |
|---|---:|---:|---:|
| `m2_l8_t30` | model2 | 8 bp | 30% |
| `m2_l10_t60` | model2 | 10 bp | 60% |

Additional presets:

| Preset | Configurations | Use |
|---|---|---|
| `genome-balanced` | `m2_l8_t30 + m2_l10_t60` | default genome-scale analysis |
| `genome-strict` | `m2_l10_t60` | shorter, high-confidence candidate lists |
| `shin-publication` | `m2_l6_t30 + m2_l10_t60` | Shin validation and publication reproduction |

Overlapping hits are merged into candidate regions. For each candidate the tool
computes:

- Z-DNA Hunter support across the two configurations
- raw Z-DNA Hunter score
- candidate length and overlap statistics
- nearest TSS distance
- nearest-TSS context score
- optional overlap annotations from user-supplied genomic or epigenomic BED tracks
- fuzzy rule activations and final fuzzy score

Three operating modes are available:

| Mode | Purpose |
|---|---|
| `balanced` | default operating point, closest to the ZDNABERT comparison |
| `moderate` | higher specificity, moderate loss of recall |
| `strict` | strongest specificity mode for shorter candidate lists |

For genome-scale runs, candidates are also filtered by `--min-score`. If not
specified, the default is `30` for genome presets. The `shin-publication` preset
uses the original publication threshold by default.

## Shin validation reference

The fuzzy layer was tuned and validated against experimentally verified Shin
loci. The operating points used by this script correspond to:

| Model | Recall | Specificity | Balanced accuracy |
|---|---:|---:|---:|
| ZDNABERT HG18 threshold 0.25 | 0.877 | 0.890 | 0.880 |
| Fuzzy detector, balanced | 0.859 | 0.904 | 0.882 |
| Fuzzy detector, moderate specificity | 0.840 | 0.945 | 0.892 |
| Fuzzy detector, strict specificity | 0.830 | 0.973 | 0.901 |

Genome-scale generation keeps the runtime profile of Z-DNA Hunter plus a small
post-processing step, while the fuzzy layer provides an adjustable decision
threshold without training a deep model for every run.

## Installation

Create an environment and install the IBP API client:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

The script does not store credentials. You can pass them on the command line or
use environment variables:

```bash
export IBP_EMAIL="your.email@example.org"
export IBP_PASSWORD="your-password"
```

## Input formats

### FASTA

The FASTA file may contain one sequence, multiple chromosomes, or a larger
genome assembly:

```text
>chr2L
ACGT...
>chr2R
ACGT...
```

To process only selected records, use `--chromosomes chr2L,chr2R`. To process a
subregion from a multi-FASTA file, use `--region chr2L:100000-200000`.

### TSS annotation

The tool accepts common SGA/BED/CSV/TSV-like files. Headered files should contain
a chromosome column and either a TSS/position column or start/end coordinates.
Headerless inputs are interpreted as one of:

```text
chr2L  12345  +
chr2L  12000  12100  gene1  0  +
```

Coordinates are treated as 0-based by default. Use `--tss-coordinate-base 1` if
single-position TSS values are 1-based.

### Optional annotation tracks

Candidate regions can be annotated against BED-like genomic context tracks and
epigenomic tracks. These tracks are exported as overlap columns and do not change
the fuzzy score unless the model is explicitly retuned.

```bash
python zdna_fuzzy_detector.py \
  --fasta data/genome.fa \
  --tss data/tss.sga \
  --annotation-bed promoters=tracks/promoters.bed \
  --annotation-bed exons=tracks/exons.bed \
  --epigenomic-bed ATAC=tracks/atac_peaks.narrowPeak \
  --output results/zdna_annotated.csv
```

## Examples

Balanced CSV output for a whole genome:

```bash
python zdna_fuzzy_detector.py \
  --fasta data/genome.fa \
  --tss data/tss.sga \
  --preset genome-balanced \
  --mode balanced \
  --min-score 30 \
  --output results/zdna_balanced.csv
```

Strict bedGraph output for one chromosome:

```bash
python zdna_fuzzy_detector.py \
  --fasta data/genome.fa \
  --tss data/tss.sga \
  --chromosomes chr2L \
  --preset genome-strict \
  --mode strict \
  --format bedgraph \
  --output results/chr2L.strict.bedgraph
```

Moderate mode on a selected interval:

```bash
python zdna_fuzzy_detector.py \
  --fasta data/genome.fa \
  --tss data/tss.sga \
  --region chr2L:1000000-1500000 \
  --preset genome-balanced \
  --mode moderate \
  --min-score 30 \
  --output results/chr2L_1Mb_1_5Mb.moderate.csv
```

Publication-like high-recall run for Shin-style validation:

```bash
python zdna_fuzzy_detector.py \
  --fasta data/genome.fa \
  --tss data/tss.sga \
  --chromosomes chr2L \
  --preset shin-publication \
  --mode balanced \
  --output results/chr2L.shin_publication_like.csv
```

Custom Z-DNA Hunter configuration:

```csv
config_id,model,min_sequence_size,threshold,score_gc,score_gtac,score_at,score_oth
custom_l8_t40,model2,8,40,2,1,0.5,0
custom_l10_t65,model2,10,65,2,1,0.5,0
```

Run with the custom configuration:

```bash
python zdna_fuzzy_detector.py \
  --fasta data/genome.fa \
  --tss data/tss.sga \
  --hunter-config configs/custom_hunter.csv \
  --mode balanced \
  --min-score 30 \
  --output results/custom_hunter_zdna.csv
```

Use unique `config_id` values for every experimental setting. When re-running
experiments with changed parameters, also use a new `--run-name` or add
`--no-reuse`.

Write all candidates, including candidates below the selected mode and
`--min-score` thresholds:

```bash
python zdna_fuzzy_detector.py \
  --fasta data/chr2L.fa \
  --tss data/tss.sga \
  --preset genome-balanced \
  --mode balanced \
  --include-all \
  --output results/chr2L.all_candidates.csv
```

## Outputs

CSV output uses 0-based half-open coordinates and contains:

- `chrom`, `start`, `end`, `candidate_id`
- `preset`, `mode`, `min_score`
- `prediction` and `fuzzy_score`
- fuzzy components: signal, context, evidence, bias, feasibility
- TSS context: `tss_distance_bp`, `tss_proximal`,
  `tss_proximity_bin`, `tss_context_label`
- Z-DNA Hunter grid features and per-configuration hit flags
- optional `annotation_*` and `epigenomic_*` overlap columns when BED tracks are supplied
- top active fuzzy rules

bedGraph output contains:

```text
chrom  start  end  fuzzy_score
```

By default only candidates accepted by the selected mode are written. Use
`--include-all` to export every merged candidate.

Each run also writes a JSON summary next to the output file:

```json
{
  "preset": "genome-balanced",
  "mode": "balanced",
  "min_score": 30.0,
  "total_candidates": 1234,
  "selected_candidates": 431,
  "mean_fuzzy_score": 22.41,
  "detection_class_counts": {
    "stable_parameter_hit": 911,
    "strict_consensus": 323
  }
}
```

## Validation and reproducibility

The standard workflow is always:

```text
FASTA + TSS -> Z-DNA Hunter API -> fuzzy detector -> CSV/bedGraph
```

For model validation, regression tests or publication reproducibility, the
script can also re-score an already prepared candidate feature table without
running the API:

```bash
python zdna_fuzzy_detector.py \
  --score-feature-table data/shin_advanced_predictions.csv \
  --preset shin-publication \
  --mode balanced \
  --include-all \
  --output results/shin_rescored.csv
```

This mode is not required for normal genome, chromosome or interval analysis.
It is useful when the candidate feature table is fixed and the goal is to
verify that the fuzzy scoring layer has not changed.

## Notes

- Credentials are read from `--email` / `--password` or from `IBP_EMAIL` /
  `IBP_PASSWORD`.
- Raw Z-DNA Hunter exports are kept under `zdna_fuzzy_runs/raw_exports` for
  reproducibility.
- Temporary FASTA slices are removed after a successful run unless
  `--keep-intermediate` is used.
- The current fuzzy model uses Z-DNA Hunter signal and nearest-TSS context. CpX
  features are not required by this tool.
- Legacy feature-table inputs with `promoter_overlap`, `regulatory_marks`,
  `repeat_overlap_pct` or `primer_uniqueness` are accepted for reproducibility,
  but newly written CSV files use neutral column names:
  `tss_proximal`, `context_support_bin`, `heuristic_bias_score` and
  `candidate_uniqueness_score`.
- `--score-feature-table` reproduces the publication scoring when the same
  candidate feature columns are supplied.
