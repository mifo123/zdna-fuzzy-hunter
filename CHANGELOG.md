# Changelog

## 1.1.0

- Added a fully local Z-DNA Hunter backend and retained the optional IBP API
  backend.
- Exposed the complete frozen fuzzy model in machine-readable form.
- Added coordinate-normalization and local-scanner regression tests.
- Added reproducible frozen-model checks for U2OS GSE290662 and rice
  GSE252954.
- Added the public 391-record HG Shin FASTA source, explicit accounting for six
  excluded non-primary contigs, locally rebuilt Hunter features, exact expected
  metrics, and an automated analysis-reproduction test.
- Corrected the source FASTA coordinate conversion from 1-based closed to
  0-based half-open; the correction changes one positive locus from FN to TP.
- Documented the exact dm6 runtime configuration: `genome-balanced` preset,
  `strict` mode and minimum score 50.
