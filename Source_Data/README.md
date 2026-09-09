# Source Data

> **Final validated status (2026-07-29):** The clustering and stability
> simulations, the PBMC RNA-only correction, Pancreas, and Bhattacherjee have
> all been updated or reconfirmed
> from completed formal reruns. The real-data PCP results use loadings estimated
> from the PCP low-rank component and the common score rule `X %*% V`. Files
> combining the three real datasets use the validated result for every dataset.

This directory contains curated CSV files underlying the current manuscript
and Supplementary Material. `figure_manifest.csv` maps each item to its
source-data file and generating script.

## Main manuscript

| Item | Source-data file(s) |
|---|---|
| Figure 1 | Manual workflow artwork; no numerical Source Data. |
| Figure 2 | `Figure02_weight_functions_source_data.csv` |
| Figure 3 | `Figure03_simulation1_gmm_pc10_clustering_ari_nmi_source_data.csv`; `simulation_figures_input_alpha005.csv` |
| Figures 4–5 | `subspace_figures_input_alpha005.csv` |
| Figure 6 | `Figure06_realdata_gmm_pc10_ari_source_data.csv`; this file retains the outer `Replicate` and inner `ClusterRep` identifiers. The two main panels plot 20 within-subset means for FULL-vs-subset and subset-vs-labels; the FULL-vs-labels mean and range from 10 clustering runs on one FULL embedding are overlaid as references in the label-agreement panel. PCP scores use `X %*% V`. |
| Figure 7 | `Figure07_realdata_subspace_pc10_source_data.csv`; unchanged by the PCP score correction because the loading vectors are unchanged. |
| Figure 8 | `Figure08_runtime_scaling_source_data.csv`; replicate-level values in `runtime_benchmark_results_matched.csv` and the unduplicated seven-setting summary in `runtime_benchmark_summary_matched.csv` |
| Figure 9 | `Figure09_method_summary_source_data.csv` (author-curated qualitative ratings); `Figure09_evidence_table.csv` (60 method-dimension rows with quantitative evidence derived from the upstream simulation, real-data, and Figure 8 runtime Source Data) |

## Supplementary Material

| Item | Source-data file(s) |
|---|---|
| Figure S1 | `FigureS01_simulation1_gmm_pc10_ari_with_k_source_data.csv`; `simulation_figures_input_alpha005.csv` |
| Figure S2 | `FigureS02_simulation1_kmeans_pc10_ari_nmi_source_data.csv`; `simulation_figures_input_alpha005.csv` |
| Figure S3 | `FigureS03_simulation1_louvain_pc10_ari_nmi_source_data.csv`; `simulation_figures_input_alpha005.csv` |
| Figure S4 | `FigureS04_simulation2_pcwise_pc10_source_data.csv`; replicate-level values in `FigureS04_simulation2_pcwise_pc10_replicate_source_data.csv`; individual PC1--PC20 values in `FigureS04_simulation2_pcwise_pc1_pc20_source_data.csv`; validation against the authoritative subspace results in `FigureS04_simulation2_pcwise_validation.csv` |
| Figure S5 | `FigureS05_simulation2_pcwise_rank_source_data.csv` |
| Figure S6 | `FigureS06_simulation1_pc10_vs_pc20_gmm_ari_source_data.csv` |
| Figure S7 | `FigureS07_realdata_pc10_vs_pc20_gmm_reproducibility_source_data.csv` |
| Figure S8 | `FigureS08_alpha_sensitivity_source_data.csv` |
| Table S1 | `simulation_figures_input_alpha005.csv`; `subspace_figures_input_alpha005.csv` |
| Table S2 | `Table_PcaGrid_QC_source_data.csv` |
| Table S3 | `Appendix_simulation1_dropout_pcp_ktau_paired_difference.csv` |
| Table S4 | `TableS04_subspace_clustering_correlations_source_data.csv` |

Raw third-party data and historical or diagnostic result versions are not part
of this curated directory. Compute-node hostnames are operational metadata and
have been removed from the public replicate-level runtime aggregate.
