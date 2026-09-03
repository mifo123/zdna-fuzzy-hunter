# Shin analysis data and provenance

The files in this directory separate public source sequences, author-derived
context variables, rebuilt Z-DNA Hunter features, and publication outputs. They
must not all be described as original Shin data.

## Source cohort and exclusions

`source/ZFS_yes.fasta` and `source/ZFS_no.fasta` are the HG Shin benchmark
FASTA files distributed through the public Z-DNABERT repository. That
repository links them to the HG Shin Google Drive folder. The sequences are
derived from the experiment reported by Shin et al.; the benchmark packaging is
from the later Z-DNABERT work.

The source FASTA files contain 391 records: 316 positive and 75 negative. A
sequence-level comparison against the UCSC NCBI36/hg18 reference confirms that
their headers use 1-based closed coordinates; each record also satisfies
`sequence length = end - start + 1`. The build script converts them to 0-based
half-open coordinates before analysis.

The reported analysis was restricted to `chr1`–`chr22`, `chrX`, and `chrY`.
Consequently, six records on non-primary `*_random` contigs were excluded:

- positive: `ZFS_88`, `ZFS_170`, `ZFS_185`, and `ZFS_265`;
- negative: `ZFS_67` and `ZFS_74`.

This leaves 385 analysed records: 312 positive and 73 negative. The complete
headers and reasons are stored in `shin_expected_metrics.json` and are checked
by the build process.

## Author-derived context variables

`shin_context_features.csv` is an author-generated intermediate from the
original analysis. It is not a file published by Shin et al. It preserves the
non-Hunter variables on which the frozen fuzzy model was fitted, including TSS
context and heuristic feasibility/bias fields. Some of those fields were
constructed by the earlier project workflow from Z-DNA Hunter overlays, CpX
Hunter overlaps, Ensembl 52 TSS distances, and deterministic heuristics.

Releasing this table makes the reported analysis exactly reproducible, but it
does not turn those engineered variables into independent experimental
measurements. This limitation must be retained in the manuscript and response
to reviewers.

## Rebuilt analysis table

`shin_analysis_features.csv` is generated, not downloaded. Run:

```bash
python validation/build_shin_analysis_features.py
```

The build verifies the 391-to-385 accounting, converts source coordinates, and
recomputes the two Z-DNA Hunter configurations with the local implementation:

- `m2_l6_t30`: model 2, minimum length 6 bp, threshold 30%;
- `m2_l10_t60`: model 2, minimum length 10 bp, threshold 60%.

The corrected local reconstruction changes one positive locus (`ZFS_104`) from
FN to TP in all three operating modes compared with the earlier API-derived
table. The expected metrics were updated rather than concealing that correction.

`supplementary/Supplementary_Table_S1_Shin_per_locus.csv` is a compact
publication output and is not a scoring input.

## Public sources

- Shin S-I et al. *DNA Research* 2016;23:477–486.
  <https://doi.org/10.1093/dnares/dsw031>
- Umerenkov D et al. *Life Science Alliance* 2023;6:e202301962.
  <https://doi.org/10.26508/lsa.202301962>
- Z-DNABERT code and HG Shin data links:
  <https://github.com/mitiau/Z-DNABERT>
- Public HG Shin folder linked by Z-DNABERT:
  <https://drive.google.com/drive/folders/1-3Ntyyjp-JfJ_V2ZXORedCDihAgckRQV>
- UCSC NCBI36/hg18 reference sequence:
  <https://hgdownload.soe.ucsc.edu/goldenPath/hg18/bigZips/>
- Ensembl release 52 TSS collection for NCBI36/hg18:
  <https://epd.expasy.org/mga/hg18/ensembl52/ensembl52.html>

The author-created Zenodo dataset is not an input to this workflow.

## SHA-256 checksums

```text
529cee5f3972c6c038c108e04088ff6ac1c4c234c3ff35ad0563801ddd9c17ad  source/ZFS_yes.fasta
c858f5875b9768e60b3c6ab4857005ac9cc844242af1061607aba5f36c437f2d  source/ZFS_no.fasta
f3f831a45f917c1c38329312e1f4ffe03724637f17dba959956cb2662e6715d0  shin_context_features.csv
94aaba1aa86cea43df2ddf09cb34103edd0d8a56a5475d87117b693d249a9ef1  shin_analysis_features.csv
```

`validation/reproduce_shin_benchmark.py` verifies these checksums, rebuilds the
analysis table in its output directory, scores all three modes, and checks the
confusion matrices and metrics in `shin_expected_metrics.json`.
