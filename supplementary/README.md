# Supplementary files

- **Table S1** gives all 385 Shin et al. loci, predictions at all three operating points, raw Z-DNA Hunter features and explicit rank shifts.
- **File S2** is a machine-readable specification of every membership function, component coefficient, rule, threshold and coordinate convention.
- **Table S3** separates quantitative comparisons from tools discussed qualitatively.
- **Table S4** reports the local runtime benchmark by stage.
- **Table S5** reports the independent GSE290662 U2OS WT chr22 validation by fuzzy-score decile.
- **Data S6** contains the deterministic 20,000-candidate sample with both Z22-IP and matched-input signals.
- **Table S7** reports the independent original-data GSE252954 rice chr1 validation by fuzzy-score decile.
- **Data S8** contains the deterministic 20,000-candidate rice sample with ZIP-seq, input, IgG and CUT&Tag signals.

Coordinates are 0-based, half-open `[start, end)`. Ranks use descending score with candidate identifier as the deterministic tie-breaker. Both independent analyses use the frozen Shin model, 2,000 candidates per decile and exact mean BigWig coverage in centred 500-bp windows. U2OS uses seed 20260803; rice uses seed 20260805 and the original DOI-associated GEO tracks.
