# Reproducibility notes

## Source of reported results

The curated manuscript Source Data are deposited with the paper rather than
duplicated in this code repository. They are the numerical basis for the
reported figures and Supplementary tables. After downloading them into
`Source_Data/`, the supplied figure scripts regenerate the corresponding
figures and derived summaries.

## Simulation settings

- Simulation 1 contains 21 tasks.
- Each Simulation 1 task has 10 independent repeats.
- Simulation 1/2 use corrected doublet generation from unmodified singlet donors.
- Splatter clean-base generation uses `dropout.type = "none"`.
- Main proposed-method calibration uses `alpha = 0.05`.
- `alpha = 0.10` and `alpha = 0.20` are used only as sensitivity analyses.
- For Simulation 1 k-means summaries, each simulation replicate contains 20
  k-means runs with `nstart = 10`; the replicate-level value is their mean.
- Louvain analyses standardize the retained PC coordinates before constructing
  the SNN graph; GMM and k-means use the retained PC scores directly.

## Real-data settings

- Datasets: PBMCs, Pancreas, and Bhattacherjee.
- PBMCs: 10x Genomics `pbmc_1k_protein_v3` (713 cells; Cell Ranger 3.0.0),
  using the official k-means six-cluster assignments as auxiliary reference
  labels rather than manually curated biological cell-type annotations.
  The raw 10x matrix is filtered to `feature_type == "Gene Expression"` before
  library-size normalization, log transformation, HVG selection, and scaling.
  In the formal input this retains 33,538 RNA features and excludes 17
  Antibody Capture features; all selected 2,000 HVGs are Gene Expression rows.
- Pancreas: GEO GSM2230759 from GSE84133 (3,605 cells), using the 14-class
  `assigned_cluster` annotation in the official processed UMI count table.
- Bhattacherjee: the official GSE124952 processed expression matrix and
  metadata (35,360 matched cells), using the eight-class `CellType` field. The
  24,822-cell derivative reported in the SCENA benchmark is not used here.
- HVGs: 2,000 for all three datasets.
- Main dimensionality: PC10.
- Sensitivity dimensionality: PC20.
- Main clustering method: GMM.
- Supplementary clustering methods: k-means and Louvain.
- The 2,000-HVG set and gene-wise scaling parameters are estimated from the FULL
  dataset and then held fixed for all subsets within that dataset.
- All real-data methods return loading vectors `V`, and all clustering and
  annotation-agreement evaluations use scores `X %*% V` from the same
  preprocessed matrix `X`. PCP estimates `V` by applying PCA to the PCP
  low-rank component `L`, but its evaluation scores are also `X %*% V`.

## Seed policy

Simulation scripts use a default master seed of `12345` unless another
`--base-seed` is supplied. Simulation seed/provenance columns such as
`MasterSeed`, `CleanSeed`, `PerturbationSeed`, `BackgroundDropoutSeed`, and
`DropoutSweepSeed` are preserved in the formal result tables.

Real-data scripts use a default `BASE_SEED = 12345`. Repeated subset analyses
derive subset and clustering seeds from the replicate index, method index, and
inner clustering replicate index. Exact formulas are implemented in the
real-data scripts.

## Computational environment

Large-scale simulations and real-data analyses were run as CPU-based jobs on
the Trubas institutional data-analysis server managed by Altair Grid Engine.
GPU acceleration was not used.

Server resources used in the formal runs:

| Analysis type | OpenMP slots per job | Memory request |
|---|---:|---:|
| Simulation 1/2 formal jobs | 1 | 8 GB per slot |
| Simulation aggregation jobs | 1 | 4--8 GB per slot |
| PBMCs real-data pairwise jobs | 5 | 8 GB per slot |
| Pancreas real-data pairwise jobs | 5 | 8 GB per slot |
| Bhattacherjee real-data pairwise jobs | 5 | 60 GB per slot |

Local figure generation and manuscript assembly were recorded under R 4.5.0 on
macOS. The recorded environments are provided in `environment/`, including the
formal-server session and package records and the figure-generation session and
package records.

## Runtime interpretation

The runtime benchmark uses gene scaling at fixed `n = 2000` and cell scaling at
fixed `d = 2000`, with common values `500, 1000, 2000, 3000`. The shared
`n = 2000, d = 2000` anchor is generated once per replicate. Five independent
3000-cell by 3000-gene master simulations are used, and all plotted inputs are
nested subsets within replicate. Each method is run in an independent
single-threaded R process. The timer includes method fitting and construction
of the 20-PC score matrix and excludes simulation, subsetting, and
preprocessing. Plotted points are medians and intervals are interquartile
ranges across the five replicate-level inputs.
The server run recorded R 4.4.1, Splatter 1.30.0, rrcov 1.7-7, rpca 0.2.3,
and RSpectra 0.16.2; these values are also retained in the replicate-level
Source Data.

Runtime values are empirical and implementation-dependent. They should be used
to interpret relative scalability under the recorded implementation and compute
environment, not as implementation-independent complexity guarantees.
