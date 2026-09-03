# Supplementary files

- **Table S1** gives the 385 primary-chromosome HG Shin analysis loci,
  predictions at all three operating points, rebuilt Z-DNA Hunter features and
  explicit rank shifts. The six excluded `*_random`-contig records are listed
  under `validation/data/`.
- **File S2** is a machine-readable specification of every membership function, component coefficient, rule, threshold and coordinate convention.
- **Table S3** separates quantitative comparisons from tools discussed qualitatively.
- **Table S4** reports the local runtime benchmark by stage. The recorded run
  used the `genome-balanced` preset (`m2_l8_t30` and `m2_l10_t60`), `strict`
  mode and minimum score 50.
- **Table S5** reports the independent GSE290662 U2OS WT chr22 validation by fuzzy-score decile.
- **Data S6** contains the deterministic 20,000-candidate sample with both Z22-IP and matched-input signals.
- **Table S7** reports the independent original-data GSE252954 rice chr1 validation by fuzzy-score decile.
- **Data S8** contains the deterministic 20,000-candidate rice sample with ZIP-seq, input, IgG and CUT&Tag signals.

Coordinates are 0-based, half-open `[start, end)`. Ranks use descending score with candidate identifier as the deterministic tie-breaker. Both independent analyses use the frozen Shin model, 2,000 candidates per decile and exact mean BigWig coverage in centred 500-bp windows. U2OS uses seed 20260803; rice uses seed 20260805 and the original DOI-associated GEO tracks.

The Shin analysis is reproduced from the versioned public FASTA source,
author-derived context table, and expected metrics under `validation/data/`;
these repository files are not part of the numbered supplementary sequence
because Supplementary Table S1 is the compact per-locus publication output.
