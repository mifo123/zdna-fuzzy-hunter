# Validation and benchmark reproduction

## Shin benchmark

The exact reported Shin benchmark can be reproduced without network access or
an API account:

```bash
python validation/reproduce_shin_benchmark.py
```

The script uses the released 385-locus feature table in `validation/data/`,
runs the public CLI in all three operating modes and verifies every confusion
matrix and reported performance metric. See `validation/data/README.md` for
input provenance, coordinate conventions and the SHA-256 checksum.

This is the canonical benchmark-reproduction command. Supplementary Table S1
is a publication output and is not accepted as a substitute input table.

## Independent experimental validation

This workflow reproduces the independent, frozen-model check reported in the
revised manuscript. It uses the human U2OS Z-DNA ChIP-seq series **GSE290662**.

## Inputs

- hg19 chromosome 22 FASTA from UCSC;
- GENCODE v19 annotation GTF;
- WT U2OS input tracks `GSM8818982` and `GSM8818983`;
- WT U2OS Z22-IP tracks `GSM8818990` and `GSM8818991`.

The four BigWig files are the GEO-processed, spike-in-normalized tracks. The
core detector remains standard-library-only; this validation helper separately
requires `pyBigWig`:

```bash
python -m pip install -r requirements-validation.txt
```

## Reproduction

Create the chr22 TSS table:

```bash
python validation/prepare_gencode_v19_tss.py \
  --gtf gencode.v19.annotation.gtf.gz \
  --chromosome chr22 \
  --output gencode_v19_chr22_tss.csv
```

Run the already-fixed Shin publication model without U2OS tuning:

```bash
python zdna_fuzzy_detector.py \
  --fasta hg19_chr22.fa \
  --tss gencode_v19_chr22_tss.csv \
  --preset shin-publication \
  --mode balanced \
  --include-all \
  --hunter-backend local \
  --output u2os_chr22_predictions.csv
```

Measure exact mean signal in fixed 500-bp windows and aggregate the deterministic
score-decile sample:

```bash
python validation/u2os_chr22_validation.py \
  --predictions u2os_chr22_predictions.csv \
  --input-rep1 GSM8818982_Input_raw_CT_1.bw \
  --input-rep2 GSM8818983_Input_raw_CT_2.bw \
  --ip-rep1 GSM8818990_IP_Z22_raw_CT_1.bw \
  --ip-rep2 GSM8818991_IP_Z22_raw_CT_2.bw \
  --output-table Supplementary_Table_S5_U2OS_validation_deciles.csv \
  --output-sample Supplementary_Data_S6_U2OS_chr22_signal_sample.csv.gz \
  --output-summary u2os_chr22_validation_summary.json
```

Sampling uses seed `20260803` and 2,000 candidates per fuzzy-score decile. An
experimental-positive window must have Z22-IP greater than matched input in
both replicates and mean Z22-IP in the top 10% of sampled windows. The criterion
was fixed before inspecting the result. It is a focused enrichment check, not a
replacement for whole-genome peak calling or a fully independent labelled
classification benchmark.

## Original-data rice validation from the reviewer-specified DOI

This second workflow uses the original GEO series **GSE252954** associated with
**doi:10.1111/pbi.14585**. Inputs are
the official MSU/UGA Release 7 chromosome FASTA and GFF3 annotation, ZIP-seq
Z22-IP replicates `GSM8010917`/`GSM8010918`, input `GSM8010919`, IgG
`GSM8010920`, and CUT&Tag replicates `GSM8010921`/`GSM8010922`.

Create the rice TSS table and run the frozen model on chromosome 1:

```bash
python validation/prepare_msu7_tss.py \
  --gff3 osa1_r7.all_models.gff3.gz \
  --chromosome Chr1 \
  --output msu7_chr1_tss.csv

python zdna_fuzzy_detector.py \
  --fasta osa1_r7.asm.chrs.fa \
  --chromosomes Chr1 \
  --tss msu7_chr1_tss.csv \
  --preset shin-publication \
  --mode balanced \
  --include-all \
  --hunter-backend local \
  --output rice_chr1_predictions.csv
```

Then measure the six original BigWig tracks:

```bash
python validation/rice_chr1_validation.py \
  --predictions rice_chr1_predictions.csv \
  --input GSM8010919_ZDNA-Input.bw \
  --igg GSM8010920_ZDNA-igG.bw \
  --ip-rep1 GSM8010917_ZDNA-IP-1.bw \
  --ip-rep2 GSM8010918_ZDNA-IP-2.bw \
  --cut-tag-rep1 GSM8010921_ZDNA-CUT_Tag-1.bw \
  --cut-tag-rep2 GSM8010922_ZDNA-CUT_Tag-2.bw \
  --output-table Supplementary_Table_S7_rice_validation_deciles.csv \
  --output-sample Supplementary_Data_S8_rice_chr1_signal_sample.csv.gz \
  --output-summary rice_chr1_validation_summary.json
```

Sampling uses seed `20260805` and 2,000 candidates per decile. A primary-positive
window requires both IP replicates to exceed both input and IgG and mean IP to
fall in the top 10% of sampled windows. CUT&Tag is a secondary continuous check.
The fixed result is 16.8% versus 2.8% in the highest and lowest deciles (risk
ratio 6.0), but all sampled highest-decile candidates are within 20 kb of a TSS.
