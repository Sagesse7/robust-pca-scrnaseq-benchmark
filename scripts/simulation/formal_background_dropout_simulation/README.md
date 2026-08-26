# Simulation workflows

This directory contains the clustering simulation (`simulation1` in the code)
and stability simulation (`simulation2`).

## Design

- Splatter base counts are generated with `dropout.type = "none"`.
- Reference, synthetic-doublet, and gene-subset mean-shift conditions share
  background dropout with `dropout_mid = 0`.
- The dropout sweep applies each specified midpoint directly to the base counts.
- The parameter grid's legacy `meanlog` field stores the median multiplicative
  factor `m`; the code uses `meanlog = log(m)` and `sdlog = 0.3`.
- Generalized K's tau cutoffs use empirical pairwise-distance quantiles. The
  primary analysis uses `alpha = 0.05`; `0.10` and `0.20` are sensitivity settings.

The default parameter grid is `config/sim_params_1_v3.xlsx`. For example:

```bash
Rscript scripts/simulation1/run_simulation1.R \
  --task-index=1 --repeats=10 --output-dir=/path/to/output/simulation1_task01
```

Use `--param-file` to supply another grid. Full runs are computationally
intensive; scheduler-specific wrappers are not included.
