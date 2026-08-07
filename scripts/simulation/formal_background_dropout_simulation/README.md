# Formal background-dropout simulation workflows

This directory provides the clustering simulation (`simulation1` in the code)
and stability simulation (`simulation2` in the code), using pairwise
empirical-quantile calibration and an explicit common background-dropout design.

Design summary:

- Base data are generated with Splatter using `dropout.type = "none"`.
- For the reference, synthetic-doublet, and gene-subset mean-shift conditions,
  a common background dropout level is applied with `dropout_mid = 0`.
- For the dropout sweep, the specified `dropout_mid` value is applied directly
  to the no-dropout base data. Dropout is not stacked twice.
- Proposed generalized K's tau variants use pairwise empirical-quantile cutoffs.
  The main analysis uses `alpha = 0.05`; `0.10` and `0.20` are sensitivity
  settings only.

The simulation scripts can be invoked directly. From this directory, the first
clustering-simulation task can be run as:

```bash
Rscript scripts/simulation1/run_simulation1.R \
  --task-index=1 --repeats=10 --output-dir=/path/to/output/simulation1_task01
```

Use `config/sim_params_1_v3.xlsx` as the default clustering-simulation parameter grid,
or supply an alternative with `--param-file`. Full runs are computationally
intensive. Scheduler-specific wrappers are intentionally omitted from the
public repository, because queue names and resource directives are local to an
individual computing system.
