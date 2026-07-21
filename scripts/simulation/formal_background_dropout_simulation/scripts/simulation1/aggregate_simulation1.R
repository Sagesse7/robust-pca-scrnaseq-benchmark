suppressPackageStartupMessages(library(data.table))

script_dir <- dirname(normalizePath(sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)[1])))
release_dir <- dirname(dirname(script_dir))
source(file.path(release_dir, "R", "simulation_utils.R"))
cli <- parse_cli()
input_root <- if (is.null(cli[["input-root"]])) file.path(release_dir, "results", "simulation1") else cli[["input-root"]]
output_dir <- if (is.null(cli[["output-dir"]])) {
  file.path(release_dir, "results", "simulation1_aggregate", format(Sys.time(), "%Y%m%d_%H%M%S"))
} else cli[["output-dir"]]
nonempty_output_guard(output_dir)

read_all <- function(filename) {
  files <- list.files(input_root, filename, recursive = TRUE, full.names = TRUE)
  if (!length(files)) stop("No ", filename, " files found.")
  rbindlist(lapply(files, fread), fill = TRUE)
}

assert_unique_keys <- function(x, keys, label) {
  missing_keys <- setdiff(keys, names(x))
  if (length(missing_keys)) {
    stop(label, " is missing key columns: ", paste(missing_keys, collapse = ", "))
  }
  duplicates <- x[, .N, by = keys][N > 1L]
  if (nrow(duplicates)) {
    preview <- paste(capture.output(print(head(duplicates, 10L))), collapse = "\n")
    stop(label, " contains duplicate result keys. Remove or explicitly select a run before aggregation.\n", preview)
  }
  x
}

result_keys <- c("ParamID", "Replicate", "Method", "Alpha")
metrics <- assert_unique_keys(read_all("simulation1_metrics.csv"), result_keys, "Simulation 1 metrics")
runtime <- assert_unique_keys(read_all("simulation1_runtime_failure.csv"), result_keys, "Simulation 1 runtime diagnostics")
metric_cols <- grep("^(ARI|NMI|k_)", names(metrics), value = TRUE)
summary <- metrics[, c(
  lapply(.SD, mean, na.rm = TRUE),
  setNames(lapply(.SD, sd, na.rm = TRUE), paste0(metric_cols, "_sd"))
), by = .(ParamID, Method, Calibration, Alpha), .SDcols = metric_cols]

data_diagnostic_cols <- intersect(c(
  "BaseID", "ParamID", "Replicate", "Condition",
  "MasterSeed", "CleanSeed", "PerturbationSeed",
  "CodeVersion", "DoubletDonorDesign",
  "NoiseDesign", "BackgroundDropoutMid", "BackgroundDropoutSeed",
  "BackgroundZeroFraction", "DropoutSweepSeed", "NoiseDesignDiagnostic",
  "NCells", "InputGenes", "AnalysisGenes",
  "DOverN", "BaseZeroFraction", "ConditionZeroFraction",
  "GroupCountMin", "GroupCountMax", "GroupCountVector",
  "AddedZeroFractionAllEntries", "AddedZeroFractionAmongNonzero",
  "ChangedEntryFraction", "ContaminatedCells", "ContaminatedCellFraction"
), names(metrics))
method_diagnostic_cols <- intersect(c(
  "BaseID", "ParamID", "Replicate", "MasterSeed", "CleanSeed",
  "PerturbationSeed", "CodeVersion", "DoubletDonorDesign",
  "NoiseDesign", "BackgroundDropoutMid",
  "BackgroundDropoutSeed", "DropoutSweepSeed",
  "Method", "Alpha", "RuntimeSec", "success",
  "failure_message", "EigenBackend", "Q1", "Q2", "Q3", "Q3star",
  "grid_scale_method", "grid_computed_k", "grid_orthogonality_error",
  "pcp_lambda", "pcp_converged", "pcp_iterations", "pcp_final_delta",
  "pcp_relative_reconstruction_error", "pcp_low_rank_estimate",
  "pcp_sparse_fraction"
), names(runtime))
data_diagnostics <- unique(metrics[, ..data_diagnostic_cols], by = c("ParamID", "Replicate"))
method_diagnostics <- runtime[, ..method_diagnostic_cols]

fwrite(metrics, file.path(output_dir, "simulation1_background_metrics_long.csv"))
fwrite(summary, file.path(output_dir, "simulation1_background_metrics_summary.csv"))
fwrite(runtime, file.path(output_dir, "simulation1_background_runtime_failure.csv"))
fwrite(data_diagnostics, file.path(output_dir, "simulation1_background_data_diagnostics.csv"))
fwrite(method_diagnostics, file.path(output_dir, "simulation1_background_method_diagnostics.csv"))
writeLines(c(
  "Simulation 1 background-dropout all-method aggregate",
  sprintf("Completed parameter settings: %d/21.", uniqueN(metrics$ParamID)),
  sprintf("Failure rows: %d.", runtime[success == FALSE, .N])
), file.path(output_dir, "SUMMARY.txt"))
