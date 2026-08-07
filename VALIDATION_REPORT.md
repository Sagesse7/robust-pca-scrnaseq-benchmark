# Code and figure validation

Validation date: 2026-08-07

This repository was validated against the manuscript, Supplementary Material,
manuscript Source Data, and final quantitative figure files. The validation
regenerated figures and derived summaries; it did not rerun the simulation or
real-data experiments.

## Manuscript-aligned analysis settings

- Main pairwise calibration: `alpha = 0.05`; `0.10` and `0.20` are sensitivity
  settings only.
- Main retained dimensionality: PC10; PC20 is a sensitivity setting.
- Real-data GMM candidate component range: `G = 1:15`.
- Real-data design: 2,000 HVGs, 20 independently sampled subsets, and 10
  clustering runs per embedding.
- Louvain shared-nearest-neighbor parameter: `k = 15`.
- PCP loadings are estimated from the low-rank component, but all downstream
  PCP scores use the common preprocessed matrix through `X %*% V`, as stated in
  the manuscript.
- PBMC preprocessing retains only features annotated as `Gene Expression`
  before normalization and HVG selection.
- Simulation and real-data method order: PCA, PcaGrid, PcaHubert, PCP, K's tau,
  Winsor, Quad, Ball, Shell, and LR.

## Checks completed

- All distributed R scripts parsed successfully.
- The formal method smoke test completed successfully for all ten methods.
- Figures 2--9 and Supplementary Figures S1--S8 were regenerated from the
  manuscript Source Data.
- All 16 regenerated quantitative PDF figures were pixel-identical to the
  corresponding manuscript figure files at 150 dpi.
- CSV files regenerated or rewritten by the figure workflow matched the
  manuscript Source Data exactly.
- Figure 1 is a manually assembled study-design schematic and is therefore not
  claimed as a code-generated figure.

## Distribution boundary

The code repository intentionally does not duplicate manuscript Source Data,
final figures, raw third-party datasets, frozen results, scheduler logs, or
institution-specific job-submission scripts. Input accessions and checksums are
documented under `data/`; the figure scripts expect the separately supplied
Source Data to be placed in `Source_Data/` before figure regeneration.
