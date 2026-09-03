# Shin benchmark data

`shin_publication_features.csv` is the immutable, 385-locus feature table used
to reproduce the reported Shin benchmark. It contains all 312 positive and 73
negative loci; no labelled locus was excluded. Coordinates are NCBI36/hg18 and
use the repository-wide 0-based, half-open convention.

The fixed paired Z-DNA Hunter configuration is:

- permissive scan: model 2, minimum length 6 bp, threshold 30%
  (`m2_l6_t30`);
- stringent scan: model 2, minimum length 10 bp, threshold 60%
  (`m2_l10_t60`).

The table includes the two scan outputs and the contextual inputs consumed by
the frozen fuzzy model. It is a released benchmark input, whereas
`supplementary/Supplementary_Table_S1_Shin_per_locus.csv` is the compact
publication output and must not be used as a scoring input.

Provenance:

- Shin S-I et al. *DNA Research* 2016;23:477–486.
  <https://doi.org/10.1093/dnares/dsw031>
- NCBI36/hg18 human reference sequence from UCSC Genome Browser.
  <https://hgdownload.soe.ucsc.edu/goldenPath/hg18/bigZips/>
- Ensembl release 52 TSS annotation (`Hum_ENSEMBL52.sga`) from the MGA/EPD
  collection. <https://epd.expasy.org/mga/hg18/ensembl52/ensembl52.html>

SHA-256:

```text
eeae4ae255adf0e5ee243305a46e17c1ba6add9fda90c120cf36b2b872e2c107  shin_publication_features.csv
```

The expected confusion matrices and metrics are recorded in
`shin_expected_metrics.json` and are checked automatically by
`validation/reproduce_shin_benchmark.py`.
