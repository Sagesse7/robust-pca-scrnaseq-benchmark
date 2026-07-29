# Formal real-data pairwise-calibration scripts

These scripts implement the PBMCs, Pancreas, and Bhattacherjee analyses.

## Formal defaults

- Ten methods: PCA, PcaGrid, PcaHubert, PCP, K's tau, Winsor, Quad, Ball, Shell, and LR.
- Pairwise empirical-quantile calibration with `alpha = 0.05`.
- 2,000 HVGs for all three datasets.
- PC10 main analysis and PC20 sensitivity analysis.
- Twenty subset repeats and ten within-replicate clustering repeats.
- GMM, k-means, and Louvain outputs; GMM is the primary manuscript analysis.
- Subset sizes: PBMCs 400, Pancreas 2,000, and Bhattacherjee 10,000.

For PBMCs, the third column of the 10x `features.tsv` file is required. Raw
matrix rows are restricted to `Gene Expression` before library-size
normalization, log transformation, variance filtering, and top-2,000-HVG
selection. The script writes `feature_filter_manifest.csv` and
`selected_hvg_features.csv` so this step can be verified directly.

The full dataset determines the 2,000-HVG feature space and gene-wise scaling
parameters. Each subset is then taken from that fixed standardized matrix so
that FULL-subset loading comparisons use the same feature coordinates.

Every method is evaluated using scores `X %*% V`, where `X` is the common
preprocessed matrix and `V` contains the method's loading vectors. PCP
estimates `V` by applying PCA to its low-rank component `L`, but does not use
the reconstruction-specific scores `L %*% V`.

## Inputs

See `../../../data/README.md` from the repository root for the expected input
layout, public accessions, and byte-verified provenance records for all three
datasets. The PBMC reference labels are the official Cell Ranger k-means
six-cluster assignments for `pbmc_1k_protein_v3`; Pancreas uses the GSM2230759
`assigned_cluster` field; and Bhattacherjee uses the eight-class `CellType`
field in the official GSE124952 metadata. Raw third-party data are not included.

## Commands

From this directory:

```bash
Rscript scripts/benchmarks/run_pbmcs_pairwise_realdata.R \
  --data-root=/path/to/data-root --output-dir=/path/to/output/pbmcs

Rscript scripts/benchmarks/run_pancreas_pairwise_realdata.R \
  --data-root=/path/to/data-root --output-dir=/path/to/output/pancreas

Rscript scripts/benchmarks/run_bhattacherjee_pairwise_realdata.R \
  --data-root=/path/to/data-root --output-dir=/path/to/output/bhattacherjee
```

All scripts accept `--methods`, `--hvg`, `--alpha`, `--repeats`,
`--cluster-repeats`, `--base-seed`, and `--workers` overrides. Scheduler
wrappers are intentionally omitted because resource directives are local to an
individual computing system.

The shared `scripts/benchmarks/pairwise_calibration_core.R` contains the
canonical K's tau and generalized K's tau implementations.
