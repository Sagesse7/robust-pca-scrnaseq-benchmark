# Reproducibility notes

## Source of reported results

The curated manuscript Source Data are supplied separately with the submission
rather than duplicated in this code repository. They are the numerical basis
for the reported figures and Supplementary tables. After placing them in
`Source_Data/`, the supplied figure scripts regenerate the corresponding
figures and derived summaries.

## Simulation settings

- The clustering simulation (`simulation1` in the code) contains 21 task IDs:
  20 displayed settings and one duplicate midpoint-0 background/reference task
  retained for provenance.
- Each simulation setting has 10 independent replicates.
- The clustering and stability simulations use corrected synthetic-doublet
  generation from unmodified singlet donors.
- For gene-subset mean shifts, the parameter grid's legacy `meanlog` field
  stores the desired median multiplicative factor `m`; the simulation function
  uses `meanlog = log(m)` and `sdlog = 0.3`. The shift-magnitude sweep uses
  `m = 1.5, 2, 3, 4, 5`, and the shifted-cell-fraction sweep fixes `m = 2`.
- Splatter clean-base generation uses `dropout.type = "none"`.
- Main proposed-method calibration uses `alpha = 0.05`.
- `alpha = 0.10` and `alpha = 0.20` are used only as sensitivity analyses.
- For clustering-simulation k-means summaries, each simulation replicate contains 20
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
- The runners also calculate k-means and Louvain outputs, but the current
  manuscript does not report these real-data clustering results.
- Reference labels are used for label-agreement analyses and to set the number
  of centers in the additional k-means calculations. They are not used for
  dimensionality reduction, GMM clustering, or Louvain clustering.
- The 2,000-HVG set and gene-wise scaling parameters are estimated from the FULL
  dataset and then held fixed for all subsets within that dataset.
- All real-data methods return loading vectors `V`. Downstream scores are
  `X %*% V` except for PcaHubert, which uses the centered scores returned by
  `rrcov`, corresponding to projection after subtraction of its fitted robust
  center. This is a constant location shift: clustering comparisons are
  translation invariant, and subspace comparisons use `V` directly. PCP
  estimates `V` using an uncentered, unscaled PCA of the PCP low-rank component
  `L` (`prcomp(L, center = FALSE, scale. = FALSE)`), but its evaluation scores
  are also `X %*% V`.

## Seed policy

Simulation scripts use a default master seed of `12345` unless another
`--base-seed` is supplied. Simulation seed/provenance columns such as
`MasterSeed`, `CleanSeed`, `PerturbationSeed`, `BackgroundDropoutSeed`, and
`DropoutSweepSeed` are preserved in the formal result tables.

Real-data scripts use a default `BASE_SEED = 12345`. Repeated subset analyses
derive subset and clustering seeds from the replicate index, method index, and
inner clustering replicate index. Exact formulas are implemented in the
real-data scripts.

To reproduce the reported results, use the complete method lists and clustering
settings specified here, because stochastic fitting can depend on the execution
sequence. The formal method order, using command-line labels, was:

- Simulations: `Grid,Hubert,PCP,PCA,Tau,Winsor,Quad,Ball,Shell,LR`.
  Clustering-simulation runs used `--clustering-seed-mode=common`.
- Real data: `PCA,Grid,Hubert,PCP,Winsor,Quad,Ball,Shell,LR,Tau`, with
  `--repeats=20 --cluster-repeats=10`.

## Computational environment

Large-scale simulations and real-data analyses were run as CPU-based jobs on
the Trubas institutional data-analysis server managed by Altair Grid Engine.
GPU acceleration was not used.

Server resources used in the formal runs:

| Analysis type | OpenMP slots per job | Memory request |
|---|---:|---:|
| Clustering- and stability-simulation formal jobs | 1 | 8 GB per slot |
| Simulation aggregation jobs | 1 | 4--8 GB per slot |
| PBMCs real-data pairwise jobs | 5 | 8 GB per slot |
| Pancreas real-data pairwise jobs | 5 | 8 GB per slot |
| Bhattacherjee real-data pairwise jobs | 5 | 60 GB per slot |

Local figure generation and manuscript assembly were recorded under R 4.5.0 on
macOS. The recorded environments are provided in `environment/`, including the
formal-server session and package records and the figure-generation session and
package records.

## Subspace--clustering complementarity

Supplementary Table S4 matches clustering-simulation outcomes to stability-simulation
PC10 subspace similarity by `ParamID` and `Method` after averaging each metric
across 10 independent replicates. The shared condition with 2% shifted cells and
a shift magnitude of 2 is displayed in both mean-shift sweeps but counted once
in this analysis. The resulting 19 unique conditions include the
`dropout.mid = 0` identity reference. Spearman correlations between subspace
similarity and ARI are reported for all 190 condition--method means and again
after excluding the 19 PcaGrid means ($n=171$). The six correlations are
descriptive; no hypothesis tests or fitted regression models are used.

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
