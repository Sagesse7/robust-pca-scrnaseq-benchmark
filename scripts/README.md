# Scripts

| Subdirectory | Purpose |
|---|---|
| `figures/` | Figure scripts using the curated files in `Source_Data/`. |
| `simulation/` | Clustering- and stability-simulation pipelines (`simulation1` and `simulation2` in the code). |
| `real_data/` | Formal PBMCs, Pancreas, and Bhattacherjee pipelines. |
| `runtime/` | Matched-design runtime benchmark and aggregation code. |
| `utils/` | Shared plotting helpers. |

`figures/build_figure09_evidence.R` derives the 10-method by 6-dimension
evidence table from the delivered simulation and real-data Source Data before
`figures/plot_main_figures_02_09.R` draws Figures 2 and 9.
`figures/build_subspace_clustering_correlations.R` derives the six
full-precision Spearman coefficients reported in Supplementary Table S4 from
the condition--method means.
`figures/plot_simulation2_pcwise_from_source_data.R` regenerates Supplementary
Figures S4 and S5 from the manuscript PC-wise Source Data.
`figures/build_dropout_paired_comparison.R` derives the paired K's tau--PCP
summary reported in Supplementary Table S3.

Scheduler-specific wrappers are intentionally excluded because queue names and
resource directives are institution-specific. Full simulation and real-data
runs require the packages listed in `../environment/` and appropriate compute
resources.
