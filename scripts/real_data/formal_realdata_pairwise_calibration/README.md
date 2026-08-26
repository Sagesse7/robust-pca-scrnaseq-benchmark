# Real-data analyses

These scripts implement the PBMC, Pancreas, and Bhattacherjee analyses.

## Settings

- Methods: PCA, PcaGrid, PcaHubert, PCP, K's tau, Winsor, Quad, Ball, Shell,
  and LR.
- Primary calibration: `alpha = 0.05`; primary dimensionality: PC10.
- Each dataset uses 2,000 highly variable genes and 20 sampled subsets, with
  10 clustering runs per embedding.
- Gaussian mixture model clustering with `G = 1:15` is primary; k-means,
  Louvain, and PC20 are supplementary analyses.
- Subset sizes are 400 for PBMC, 2,000 for Pancreas, and 10,000 for
  Bhattacherjee.

Reference labels are used only for label-agreement analyses and to set the
k-means center count. PBMC features are restricted to `Gene Expression`
before preprocessing. The full dataset defines the highly variable genes and
scaling parameters used for its subsets.

Cell scores are computed as `X %*% V`, except that PcaHubert uses the centered
scores returned by `rrcov`. PCP estimates `V` from its low-rank component
with `prcomp(L, center = FALSE, scale. = FALSE)`; its evaluation scores remain
`X %*% V`. The canonical K's tau implementations are in
`scripts/benchmarks/pairwise_calibration_core.R`.

## Inputs and commands

See `../../../data/README.md` for data sources and local file layout. From
this directory, run:

```bash
Rscript scripts/benchmarks/run_pbmcs_pairwise_realdata.R \
  --data-root=/path/to/data-root --output-dir=/path/to/output/pbmcs

Rscript scripts/benchmarks/run_pancreas_pairwise_realdata.R \
  --data-root=/path/to/data-root --output-dir=/path/to/output/pancreas

Rscript scripts/benchmarks/run_bhattacherjee_pairwise_realdata.R \
  --data-root=/path/to/data-root --output-dir=/path/to/output/bhattacherjee
```

The runners accept overrides for methods, HVGs, calibration, repetitions,
seeds, and worker counts. Scheduler-specific wrappers are not included.
