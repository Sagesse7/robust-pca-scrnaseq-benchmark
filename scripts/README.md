# Scripts

| Subdirectory | Purpose |
|---|---|
| `figures/` | Figure scripts using the curated files in `Source_Data/`. |
| `simulation/` | Formal Simulation 1 and Simulation 2 pipelines. |
| `real_data/` | Formal PBMCs, Pancreas, and Bhattacherjee pipelines. |
| `runtime/` | Matched-design runtime benchmark and aggregation code. |
| `utils/` | Shared plotting helpers. |

`figures/build_figure09_evidence.R` derives the 10-method by 5-dimension
evidence table from the delivered simulation and real-data Source Data before
`figures/plot_main_figures_02_09.R` draws Figures 2 and 9.

Scheduler-specific wrappers are intentionally excluded because queue names and
resource directives are institution-specific. Full simulation and real-data
runs require the packages listed in `../environment/` and appropriate compute
resources.
