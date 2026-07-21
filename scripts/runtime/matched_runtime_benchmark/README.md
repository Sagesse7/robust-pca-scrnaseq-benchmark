# Matched-design runtime benchmark

This directory contains the formal code used for the manuscript runtime
benchmark. The benchmark uses seven unique settings formed from two matched
sweeps:

- gene scaling: `n = 2000`, `d = 500, 1000, 2000, 3000`;
- cell scaling: `d = 2000`, `n = 500, 1000, 2000, 3000`.

The `n = 2000, d = 2000` anchor is shared. Five independent replicate-level
master simulations are generated, and all settings within a replicate are
nested subsets of the same 3000-cell by 3000-gene master matrix. Each method
call uses one thread. The timer includes dimensionality-reduction fitting and
construction of the 20-PC score matrix, but excludes simulation, subsetting,
and preprocessing.

## Files

- `config/runtime_grid.csv`: seven unique settings;
- `scripts/run_runtime_task.R`: one setting-replicate-method call;
- `scripts/emit_runtime_schedule.R`: deterministic randomized call order;
- `scripts/aggregate_runtime_results.R`: validates all 350 task rows and
  reports medians and interquartile ranges;
- `scripts/plot_runtime_results.R`: converts a newly aggregated formal run to
  the matched two-panel figure and plotted summary;
- `R/figure_style.R`: runtime figure labels and plotting order;
- `../../simulation/formal_background_dropout_simulation/R/`: shared method
  and simulation implementations used by both simulation and runtime analyses.

The public consolidated replicate-level and summary results are in
`../../../Source_Data/runtime_benchmark_results_matched.csv` and
`../../../Source_Data/runtime_benchmark_summary_matched.csv`. Institution-specific
scheduler wrappers and individual job files are intentionally excluded.
