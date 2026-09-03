# HG Shin source FASTA

These two files are the HG Shin benchmark FASTA files distributed through the
public Z-DNABERT project:

- <https://github.com/mitiau/Z-DNABERT>
- <https://drive.google.com/drive/folders/1-3Ntyyjp-JfJ_V2ZXORedCDihAgckRQV>

They contain 316 positive and 75 negative records derived from the Shin et al.
Z-DNA ChIP-seq study (<https://doi.org/10.1093/dnares/dsw031>). Z-DNABERT is
described at <https://doi.org/10.26508/lsa.202301962>.

A sequence-level comparison against the UCSC NCBI36/hg18 reference confirms
that the FASTA headers encode 1-based closed intervals. Files in this repository
are line-ending-normalized to LF; their checksums are recorded in
`../shin_expected_metrics.json`.
