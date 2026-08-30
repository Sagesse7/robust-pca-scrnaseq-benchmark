# Real-data analyses

These scripts implement the PBMC, Pancreas, and Bhattacherjee analyses.

## Settings

- Methods: PCA, PcaGrid, PcaHubert, PCP, K's tau, Winsor, Quad, Ball, Shell,
  and LR.
- Primary calibration: `alpha = 0.05`; primary dimensionality: PC10.
- Each dataset uses 2,000 highly variable genes and 20 sampled subsets, with
  10 clustering runs per embedding.
- The manuscript reports Gaussian mixture model clustering with `G = 1:15`
  and PC20 sensitivity. The runners also calculate k-means and Louvain results,
  which are not reported as real-data analyses in the current manuscript.
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

See [`data/README.md`](../../../data/README.md) for input preparation. Run this
Bash block from the repository root, after installing the recorded packages.
Choose fresh output directories and retain the complete method list and
clustering settings; see [`Seed policy`](../../../docs/reproducibility_notes.md#seed-policy).

```bash
set -euo pipefail
REAL=scripts/real_data/formal_realdata_pairwise_calibration/scripts/benchmarks
DATA_ROOT=/path/to/data-root
OUT="$PWD/results/real_data"
common=(
  --data-root="$DATA_ROOT" --hvg=2000 --alpha=0.05
  --base-seed=12345 --repeats=20 --cluster-repeats=10 --workers=5
  --methods=PCA,Grid,Hubert,PCP,Winsor,Quad,Ball,Shell,LR,Tau
)
Rscript "$REAL/run_pbmcs_pairwise_realdata.R" "${common[@]}" \
  --subset-n=400 --output-dir="$OUT/pbmcs"
Rscript "$REAL/run_pancreas_pairwise_realdata.R" "${common[@]}" \
  --subset-n=2000 --output-dir="$OUT/pancreas"
Rscript "$REAL/run_bhattacherjee_pairwise_realdata.R" "${common[@]}" \
  --subset-n=10000 --output-dir="$OUT/bhattacherjee"
```

Each dataset directory contains `real_full_results_repeated_clustering.csv`,
`real_sampling_results_repeated_clustering.csv`, their summaries, and run
metadata. Retain the replicate-level CSVs: the Source Data preparation step
uses them to preserve the subset/clustering hierarchy. Commands for Figures
6, 7, and S7 are in [`scripts/README.md`](../../README.md).

These are computationally intensive jobs, especially for the FULL
Bhattacherjee data. Scheduler-specific wrappers are not included.
