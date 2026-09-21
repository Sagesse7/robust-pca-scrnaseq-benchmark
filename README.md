# Benchmarking and Generalizing Robust Principal Component Analysis for Single-Cell RNA Sequencing

This repository is the reproducibility companion for the manuscript of the
same title. It contains the analysis and figure-generation code, curated Source
Data underlying the reported figures and tables, final figure files, and the
recorded R environment.

## v1.0.2 presentation update

This version aligns the figure files and plotting scripts with the revised
manuscript. Changes concern font sizes, panel layout, axis-label orientation
and spacing, and colorbar-title alignment. Numerical Source Data, analysis
methods, and scientific results are unchanged. Figures 1 and 9 retain their
previous appearance. See [`CHANGELOG.md`](CHANGELOG.md) for details.


## Contents

| Path | Contents |
|---|---|
| `scripts/simulation/` | Clustering- and stability-simulation workflows. |
| `scripts/real_data/` | PBMC, Pancreas, and Bhattacherjee analyses. |
| `scripts/runtime/` | Matched runtime benchmark. |
| `scripts/figures/` | Figure-generation scripts. |
| `Source_Data/` | Curated CSV files underlying manuscript and Supplementary results. |
| `figures/` | Final manuscript and Supplementary figure files. |
| `data/` | Data sources, expected file layout, and checksums. |
| `environment/` | Recorded R environments. |
| `figure_manifest.csv` | Mapping from each figure or table to its Source Data and generating script. |

Raw third-party data, frozen intermediate results, server logs, and
computing-system-specific job scripts are not included.

Software-version records are in [`environment/`](environment/). The commands
below assume the required R packages have been installed.

## Analysis overview

The benchmark compares PCA, PcaGrid, PcaHubert, PCP, K's tau, Winsor, Quad,
Ball, Shell, and LR. The primary analyses use `alpha = 0.05`, the first 10
principal components, and Gaussian mixture model clustering with candidate
component counts `G = 1:15`. Additional settings are documented in
[`docs/reproducibility_notes.md`](docs/reproducibility_notes.md).

## Running the analyses

Data accessions, download links, checksums, and the required local layout are
given in [`data/README.md`](data/README.md). Analysis commands are provided in:

- [`scripts/simulation/formal_background_dropout_simulation/README.md`](scripts/simulation/formal_background_dropout_simulation/README.md)
- [`scripts/real_data/formal_realdata_pairwise_calibration/README.md`](scripts/real_data/formal_realdata_pairwise_calibration/README.md)
- [`scripts/runtime/matched_runtime_benchmark/README.md`](scripts/runtime/matched_runtime_benchmark/README.md)

To regenerate manuscript figures from the supplied CSV files in
`Source_Data/`, use the scripts described in
[`scripts/README.md`](scripts/README.md). Figures are written to `figures/`,
replacing the supplied copies. Figure 1 is a manually assembled schematic;
Figures 2--9 and Supplementary Figures S1--S8 are script-generated.

The same guide documents conversion of full-rerun outputs into Source Data.

## Validation

The release checks are summarized in
[`VALIDATION_REPORT.md`](VALIDATION_REPORT.md).

## License

Code and documentation are distributed under the BSD 3-Clause License. See
`LICENSE`. Third-party data and software retain their original terms.
