# Formal background-dropout simulation bundle

This directory provides the Simulation 1 and Simulation 2 workflows using
pairwise empirical-quantile calibration and an explicit common background
dropout design.

Design summary:

- Base data are generated with Splatter using `dropout.type = "none"`.
- For clean, doublet, and shift conditions, a common background dropout level is
  applied with `dropout_mid = 0`.
- For the dropout sweep, the specified `dropout_mid` value is applied directly
  to the no-dropout base data. Dropout is not stacked twice.
- Proposed generalized Kendall methods use pairwise empirical-quantile cutoffs
  with alpha values `0.05, 0.10, 0.20`.

The simulation scripts can be invoked directly. From this directory, the first
Simulation 1 task can be run as:

```bash
Rscript scripts/simulation1/run_simulation1.R \
  --task-index=1 --repeats=10 --output-dir=/path/to/output/simulation1_task01
```

Use `config/sim_params_1_v3.xlsx` as the default Simulation 1 parameter grid,
or supply an alternative with `--param-file`. Full runs are computationally
intensive. Scheduler-specific wrappers are intentionally omitted from the
public repository, because queue names and resource directives are local to an
individual computing system.
