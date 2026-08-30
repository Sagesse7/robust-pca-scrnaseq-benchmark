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

## Run the analyses

Run the following blocks in the same Bash session from the repository root,
using the packages recorded in `environment/`. Full runs are computationally
intensive. Use a new output location; task and aggregation scripts reject
nonempty output directories. Scheduler wrappers are not included.

```bash
set -euo pipefail
SIM=scripts/simulation/formal_background_dropout_simulation
OUT="$PWD/results"
mkdir -p "$OUT"
common=(
  --param-file="$SIM/config/sim_params_1_v3.xlsx"
  --methods=Grid,Hubert,PCP,PCA,Tau,Winsor,Quad,Ball,Shell,LR
  --repeats=10 --base-seed=12345 --kpc=20
  --alphas=0.05,0.10,0.20 --de-fac-loc=0.1 --background-dropout-mid=0
)
Rscript "$SIM/scripts/reproducibility/write_session_info.R" "$OUT/sessionInfo.txt"
```

### Clustering simulation

Run all 21 task IDs, then aggregate. Task 1 is the duplicate midpoint-0
background/reference output retained for provenance.

```bash
for task in {1..21}; do
  tag=$(printf '%02d' "$task")
  Rscript "$SIM/scripts/simulation1/run_simulation1.R" "${common[@]}" \
    --task-index="$task" --clustering-seed-mode=common --gmm-itmax=1000 \
    --output-dir="$OUT/formal/simulation1/task_$tag"
done
Rscript "$SIM/scripts/simulation1/aggregate_simulation1.R" \
  --input-root="$OUT/formal/simulation1" --output-dir="$OUT/simulation1"
```

### Stability simulation

Generate the shared reference with task 1, run perturbed tasks 2--21, then
aggregate the paired loadings. The code's `clean-reference` name denotes the
midpoint-0 background reference, not a dataset without dropout. Keep the RDS
files in place until aggregation and PC-wise analysis are complete.

```bash
Rscript "$SIM/scripts/simulation2/run_simulation2.R" "${common[@]}" \
  --mode=clean-reference --task-index=1 \
  --output-dir="$OUT/formal/simulation2_clean_reference/task_01"
for task in {2..21}; do
  tag=$(printf '%02d' "$task")
  Rscript "$SIM/scripts/simulation2/run_simulation2.R" "${common[@]}" \
    --mode=perturbed --task-index="$task" \
    --output-dir="$OUT/formal/simulation2_perturbed/task_$tag"
done
Rscript "$SIM/scripts/simulation2/aggregate_simulation2.R" \
  --clean-root="$OUT/formal/simulation2_clean_reference" \
  --perturbed-root="$OUT/formal/simulation2_perturbed" \
  --output-dir="$OUT/simulation2"
```

The aggregates include `simulation1_background_metrics_long.csv`,
`simulation1_background_metrics_summary.csv`,
`simulation2_background_similarity_long.csv`, and
`simulation2_background_similarity_summary.csv`. Check the accompanying
`SUMMARY.txt` and failure tables before export. Source Data preparation and
PC-wise commands are in [`scripts/README.md`](../../README.md).

For a lightweight implementation check before a full run:

```bash
Rscript "$SIM/tests/smoke_test_methods.R"
```

This is a limited smoke test, not a complete benchmark validation. Preserve
the formal method order and settings when reproducing reported values; see
[`Seed policy`](../../../docs/reproducibility_notes.md#seed-policy).
