suppressPackageStartupMessages(library(data.table))

script_dir <- dirname(normalizePath(sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)[1])))
release_dir <- dirname(dirname(script_dir))
source(file.path(release_dir, "R", "simulation_utils.R"))
cli <- parse_cli()
clean_root <- if (is.null(cli[["clean-root"]])) {
  file.path(release_dir, "results", "simulation2_clean_reference")
} else cli[["clean-root"]]
perturbed_root <- if (is.null(cli[["perturbed-root"]])) {
  file.path(release_dir, "results", "simulation2_perturbed")
} else cli[["perturbed-root"]]
output_dir <- if (is.null(cli[["output-dir"]])) {
  file.path(release_dir, "results", "simulation2_aggregate", format(Sys.time(), "%Y%m%d_%H%M%S"))
} else cli[["output-dir"]]
nonempty_output_guard(output_dir)

read_manifests <- function(root, filename) {
  files <- list.files(root, filename, recursive = TRUE, full.names = TRUE)
  if (!length(files)) stop("No ", filename, " files found under ", root)
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

clean_keys <- c("BaseID", "Replicate", "Method", "Alpha")
perturbed_keys <- c("ParamID", "Replicate", "Method", "Alpha")
clean <- assert_unique_keys(
  read_manifests(clean_root, "simulation2_manifest.csv"), clean_keys,
  "Simulation 2 clean-reference manifest"
)
perturbed <- assert_unique_keys(
  read_manifests(perturbed_root, "simulation2_manifest.csv"), perturbed_keys,
  "Simulation 2 perturbed manifest"
)
if (any(clean$Condition != "clean")) stop("Clean-reference root contains non-clean rows.")
if (any(perturbed$Condition == "clean")) stop("Perturbed root contains clean rows.")

runtime <- rbindlist(list(
  read_manifests(clean_root, "simulation2_runtime_failure.csv"),
  read_manifests(perturbed_root, "simulation2_runtime_failure.csv")
), fill = TRUE)
runtime <- assert_unique_keys(
  runtime, c("RunMode", "ParamID", "Replicate", "Method", "Alpha"),
  "Simulation 2 runtime diagnostics"
)

subspace_similarity <- function(A, B, k = 10L) {
  A <- qr.Q(qr(as.matrix(A)[, seq_len(min(k, ncol(A))), drop = FALSE]))
  B <- qr.Q(qr(as.matrix(B)[, seq_len(min(k, ncol(B))), drop = FALSE]))
  s <- svd(crossprod(A, B), nu = 0, nv = 0)$d
  mean(pmin(pmax(s, 0), 1))
}

clean_keys <- clean[, .(
  BaseID, Replicate, Method, Alpha,
  CleanParamID = ParamID, CleanSeedReference = CleanSeed,
  CleanGeneHash = GeneHash, CleanPath = RDS_Path
)]
pairs <- merge(
  perturbed, clean_keys,
  by = c("BaseID", "Replicate", "Method", "Alpha"),
  all.x = TRUE
)
if (anyNA(pairs$CleanPath)) {
  stop("At least one perturbed result has no matching clean reference.")
}
if (any(pairs$CleanSeed != pairs$CleanSeedReference)) {
  stop("Clean-seed mismatch between perturbed and clean-reference results.")
}

pairs[, SimilarityPC10 := mapply(function(clean_path, perturbed_path) {
  clean_obj <- readRDS(clean_path)
  perturbed_obj <- readRDS(perturbed_path)
  subspace_similarity(clean_obj$loadings_full, perturbed_obj$loadings_full, 10L)
}, CleanPath, RDS_Path)]
pairs[, SimilarityPC20 := mapply(function(clean_path, perturbed_path) {
  clean_obj <- readRDS(clean_path)
  perturbed_obj <- readRDS(perturbed_path)
  subspace_similarity(clean_obj$loadings_full, perturbed_obj$loadings_full, 20L)
}, CleanPath, RDS_Path)]
summary <- pairs[, .(
  SimilarityPC10_mean = mean(SimilarityPC10, na.rm = TRUE),
  SimilarityPC10_sd = sd(SimilarityPC10, na.rm = TRUE),
  SimilarityPC20_mean = mean(SimilarityPC20, na.rm = TRUE),
  SimilarityPC20_sd = sd(SimilarityPC20, na.rm = TRUE)
), by = .(ParamID, Condition, Method, Calibration, Alpha)]

manifest <- rbindlist(list(clean, perturbed), fill = TRUE)
fwrite(manifest, file.path(output_dir, "simulation2_background_manifest.csv"))
fwrite(pairs, file.path(output_dir, "simulation2_background_similarity_long.csv"))
fwrite(summary, file.path(output_dir, "simulation2_background_similarity_summary.csv"))
fwrite(runtime, file.path(output_dir, "simulation2_background_runtime_failure.csv"))
data_diagnostic_cols <- intersect(c(
  "RunMode", "BaseID", "ParamID", "Replicate", "Condition",
  "MasterSeed", "CleanSeed", "PerturbationSeed",
  "CodeVersion", "DoubletDonorDesign",
  "NoiseDesign", "BackgroundDropoutMid", "BackgroundDropoutSeed",
  "BackgroundZeroFraction", "DropoutSweepSeed", "NoiseDesignDiagnostic",
  "NCells", "InputGenes", "AnalysisGenes", "DOverN",
  "GroupCountMin", "GroupCountMax", "GroupCountVector",
  "BaseZeroFraction", "ConditionZeroFraction",
  "AddedZeroFractionAllEntries", "AddedZeroFractionAmongNonzero",
  "ChangedEntryFraction", "ContaminatedCells", "ContaminatedCellFraction"
), names(manifest))
method_diagnostic_cols <- intersect(c(
  "RunMode", "BaseID", "ParamID", "Replicate", "Condition",
  "MasterSeed", "CleanSeed", "PerturbationSeed",
  "CodeVersion", "DoubletDonorDesign",
  "NoiseDesign", "BackgroundDropoutMid", "BackgroundDropoutSeed",
  "DropoutSweepSeed",
  "Method", "Alpha", "RuntimeSec", "success", "failure_message",
  "EigenBackend", "Q1", "Q2", "Q3", "Q3star",
  "grid_scale_method", "grid_computed_k", "grid_orthogonality_error",
  "pcp_lambda", "pcp_converged", "pcp_iterations", "pcp_final_delta",
  "pcp_relative_reconstruction_error", "pcp_low_rank_estimate",
  "pcp_sparse_fraction"
), names(runtime))
data_diagnostics <- unique(
  manifest[, ..data_diagnostic_cols],
  by = c("RunMode", "ParamID", "Replicate", "Condition")
)
method_diagnostics <- runtime[, ..method_diagnostic_cols]
fwrite(data_diagnostics, file.path(output_dir, "simulation2_background_data_diagnostics.csv"))
fwrite(method_diagnostics, file.path(output_dir, "simulation2_background_method_diagnostics.csv"))
writeLines(c(
  "Simulation 2 background-dropout all-method aggregate",
  sprintf("Completed perturbed parameter settings: %d/20.", uniqueN(perturbed$ParamID)),
  sprintf("Clean references: %d BaseID x replicate combinations.", uniqueN(clean[, .(BaseID, Replicate)])),
  sprintf("Failure rows: %d.", runtime[success == FALSE, .N])
), file.path(output_dir, "SUMMARY.txt"))
