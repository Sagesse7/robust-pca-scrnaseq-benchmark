# Matched runtime benchmark

The benchmark uses two sweeps:

- genes: `n = 2000`, `d = 500, 1000, 2000, 3000`;
- cells: `d = 2000`, `n = 500, 1000, 2000, 3000`.

The shared `n = 2000, d = 2000` setting is generated once per replicate.
Inputs are nested subsets of five independent 3000-cell by 3000-gene matrices.
Each method uses one thread. Timing includes method fitting and 20-PC score
construction, but excludes simulation, subsetting, and preprocessing.

## Files

| Path | Purpose |
|---|---|
| `config/runtime_grid.csv` | Seven unique benchmark settings. |
| `scripts/run_runtime_task.R` | Runs one setting, replicate, and method. |
| `scripts/emit_runtime_schedule.R` | Creates the randomized call order. |
| `scripts/aggregate_runtime_results.R` | Validates and summarizes task results. |
| `scripts/plot_runtime_results.R` | Produces the runtime figure and summary. |

Scheduler-specific wrappers and job files are not included.
