# Changelog

## 1.1.0

- Added a fully local Z-DNA Hunter backend and retained the optional IBP API
  backend.
- Exposed the complete frozen fuzzy model in machine-readable form.
- Added coordinate-normalization and local-scanner regression tests.
- Added reproducible frozen-model checks for U2OS GSE290662 and rice
  GSE252954.
- Added the released 385-locus Shin benchmark input, exact expected metrics and
  an automated end-to-end reproduction test.
- Documented the exact dm6 runtime configuration: `genome-balanced` preset,
  `strict` mode and minimum score 50.
