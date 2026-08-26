# Scripts

| Subdirectory | Purpose |
|---|---|
| `simulation/` | Clustering and stability simulations. |
| `real_data/` | PBMC, Pancreas, and Bhattacherjee analyses. |
| `runtime/` | Matched runtime benchmark. |
| `figures/` | Manuscript figure and derived-summary generation. |
| `utils/` | Shared plotting helpers. |

After placing the manuscript CSV files in `Source_Data/`, run the required
figure scripts from the repository root. The main entry points are:

```bash
Rscript scripts/figures/plot_simulation1_clustering_combined.R
Rscript scripts/figures/plot_simulation_support_from_source_data.R
Rscript scripts/figures/plot_simulation2_pcwise_from_source_data.R
Rscript scripts/figures/plot_realdata_from_source_data.R
Rscript scripts/figures/plot_runtime_from_source_data.R
Rscript scripts/figures/build_figure09_evidence.R
Rscript scripts/figures/plot_main_figures_02_09.R
Rscript scripts/figures/plot_alpha_sensitivity_from_source_data.R
Rscript scripts/figures/build_dropout_paired_comparison.R
Rscript scripts/figures/build_subspace_clustering_correlations.R
```

Analysis-specific commands and settings are documented in the README within
each analysis directory. Scheduler wrappers are not included.
