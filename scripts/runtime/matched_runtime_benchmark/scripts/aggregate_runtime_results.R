#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

script_file <- sub(
  "^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
)
script_dir <- dirname(normalizePath(script_file))
bundle_root <- dirname(script_dir)

parse_cli <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    out[[parts[1]]] <- paste(parts[-1], collapse = "=")
  }
  out
}

cli <- parse_cli(commandArgs(trailingOnly = TRUE))
input_dir <- if (is.null(cli[["input-dir"]])) {
  file.path(bundle_root, "results", "tasks")
} else cli[["input-dir"]]
output_dir <- if (is.null(cli[["output-dir"]])) {
  file.path(bundle_root, "results", "aggregated")
} else cli[["output-dir"]]
grid_file <- if (is.null(cli[["grid-file"]])) {
  file.path(bundle_root, "config", "runtime_grid.csv")
} else cli[["grid-file"]]
repeats <- if (is.null(cli$repeats)) 5L else as.integer(cli$repeats)
task_id_repeats <- if (is.null(cli[["task-id-repeats"]])) {
  5L
} else {
  as.integer(cli[["task-id-repeats"]])
}

method_order <- c(
  "PCA", "Grid", "Hubert", "PCP", "Tau",
  "Winsor", "Quad", "Ball", "Shell", "LR"
)

if (!dir.exists(input_dir)) stop("Runtime task directory not found: ", input_dir)
if (!file.exists(grid_file)) stop("Runtime grid not found: ", grid_file)
grid <- fread(grid_file)

if (repeats < 1L || task_id_repeats < repeats) {
  stop("task-id-repeats must be at least as large as repeats.")
}
expected_all <- as.data.table(expand.grid(
  SettingID = grid$SettingID,
  Replicate = seq_len(task_id_repeats),
  Method = method_order,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
))
expected_all[, TaskID := .I]
expected <- expected_all[Replicate <= repeats]

files <- sort(list.files(
  input_dir, pattern = "^runtime_task_[0-9]{4}[.]csv$", full.names = TRUE
))
if (!length(files)) stop("No runtime task CSV files were found in: ", input_dir)

rows <- rbindlist(lapply(files, fread), fill = TRUE, use.names = TRUE)
required <- c(
  "TaskID", "SettingID", "Replicate", "Method", "RuntimeSec", "Success",
  "Calibration", "InputDesign", "MasterSeed", "CellSubsetSignature",
  "GeneSubsetSignature", "InputSumSquares", "FitSec", "ScoreSec", "CodeVersion",
  "MethodFileMD5", "SimulationFileMD5"
)
missing_columns <- setdiff(required, names(rows))
if (length(missing_columns)) {
  stop("Task results are missing columns: ", paste(missing_columns, collapse = ", "))
}

key <- c("TaskID", "SettingID", "Replicate", "Method")
duplicates <- rows[duplicated(rows, by = key) | duplicated(rows, by = key, fromLast = TRUE)]
if (nrow(duplicates)) {
  stop("Duplicate runtime task keys detected. Refusing silent deduplication.")
}

missing_tasks <- expected[!rows, on = key]
unexpected_tasks <- rows[!expected, on = key]
if (nrow(missing_tasks)) {
  stop("Missing ", nrow(missing_tasks), " of ", nrow(expected), " expected runtime tasks.")
}
if (nrow(unexpected_tasks)) {
  stop("Found ", nrow(unexpected_tasks), " unexpected runtime tasks.")
}
if (nrow(rows) != nrow(expected)) {
  stop("Unexpected aggregate row count after key validation.")
}
if (any(!rows$Success)) {
  failed <- rows[!Success, .(TaskID, SettingID, Replicate, Method, FailureMessage)]
  print(failed)
  stop("One or more runtime tasks failed; formal aggregation was not written.")
}
if (any(!is.finite(rows$RuntimeSec) | rows$RuntimeSec <= 0)) {
  stop("RuntimeSec must be finite and positive for every task.")
}
if (any(!is.finite(rows$FitSec) | !is.finite(rows$ScoreSec) |
        rows$FitSec <= 0 | rows$ScoreSec < 0)) {
  stop("FitSec and ScoreSec must be finite for every successful task.")
}
if (any(abs(rows$RuntimeSec - rows$FitSec - rows$ScoreSec) > 1e-6)) {
  stop("RuntimeSec must equal FitSec + ScoreSec for every task.")
}
if (uniqueN(rows$CodeVersion) != 1L ||
    uniqueN(rows$MethodFileMD5) != 1L ||
    uniqueN(rows$SimulationFileMD5) != 1L) {
  stop("Mixed code versions or source-file checksums detected.")
}

proposed <- c("Winsor", "Quad", "Ball", "Shell", "LR")
if (rows[Method %in% proposed, !all(Calibration == "pairwise_empirical_quantile")]) {
  stop("A proposed method row does not use pairwise empirical-quantile calibration.")
}
if (rows[!Method %in% proposed, !all(Calibration == "none")]) {
  stop("A reference method row has an unexpected calibration label.")
}
if (uniqueN(rows$InputDesign) != 1L ||
    unique(rows$InputDesign) != "nested_subsets_from_common_3000x3000_splatter_background_mid0") {
  stop("Runtime inputs do not have the required matched nested-subset design.")
}

# Each setting/method in a replicate must reconstruct the identical input.
input_check <- rows[, .(
  MasterSeeds = uniqueN(MasterSeed),
  InputSums = uniqueN(InputSumSquares),
  CellSignatures = uniqueN(CellSubsetSignature),
  GeneSignatures = uniqueN(GeneSubsetSignature)
), by = .(Replicate, SettingID)]
if (any(input_check$MasterSeeds != 1L | input_check$InputSums != 1L |
        input_check$CellSignatures != 1L | input_check$GeneSignatures != 1L)) {
  stop("Methods within a setting/replicate did not receive the same reconstructed input.")
}

# The n=2000 and d=2000 anchors are intentionally shared between panels.
anchor_check <- rows[Method == "PCA", .(
  cell_shared = uniqueN(CellSubsetSignature[SettingID %in% 1:4]) == 1L,
  gene_shared = uniqueN(GeneSubsetSignature[SettingID %in% c(3L, 5L, 6L, 7L)]) == 1L
), by = Replicate]
if (!all(anchor_check$cell_shared) || !all(anchor_check$gene_shared)) {
  stop("The required nested cell/gene subsets or shared anchor were not reconstructed.")
}

generalized <- c("Tau", "Winsor", "Quad", "Ball", "Shell", "LR")
generalized_check <- rows[Method %in% generalized]
if (any(!is.finite(generalized_check$GeometrySec) |
        !is.finite(generalized_check$WeightSec) |
        !is.finite(generalized_check$ScatterSec) |
        !is.finite(generalized_check$EigenSec))) {
  stop("Generalized-method component timings are incomplete.")
}

rows[, MethodOrder := match(Method, method_order)]
setorder(rows, SettingID, Replicate, MethodOrder)
rows[, MethodOrder := NULL]
summary <- rows[, .(
  N = .N,
  AnalysisGenesMin = min(AnalysisGenes),
  AnalysisGenesMedian = median(AnalysisGenes),
  AnalysisGenesMax = max(AnalysisGenes),
  MeanRuntimeSec = mean(RuntimeSec),
  SDRuntimeSec = sd(RuntimeSec),
  SERuntimeSec = sd(RuntimeSec) / sqrt(.N),
  MedianRuntimeSec = median(RuntimeSec),
  Q25RuntimeSec = quantile(RuntimeSec, 0.25, names = FALSE, type = 7),
  Q75RuntimeSec = quantile(RuntimeSec, 0.75, names = FALSE, type = 7),
  MinRuntimeSec = min(RuntimeSec),
  MaxRuntimeSec = max(RuntimeSec),
  MedianFitSec = median(FitSec),
  MedianScoreSec = median(ScoreSec),
  MedianGeometrySec = if (all(is.na(GeometrySec))) NA_real_ else median(GeometrySec, na.rm = TRUE),
  MedianWeightSec = if (all(is.na(WeightSec))) NA_real_ else median(WeightSec, na.rm = TRUE),
  MedianScatterSec = if (all(is.na(ScatterSec))) NA_real_ else median(ScatterSec, na.rm = TRUE),
  MedianEigenSec = if (all(is.na(EigenSec))) NA_real_ else median(EigenSec, na.rm = TRUE)
), by = .(
  SettingID, Panel, Cell, Gene, CellType, de_prob, DeFacLoc,
  BackgroundDropoutMid, InputDesign, MasterCell, MasterGene, Method, MethodLabel, Calibration, Alpha,
  PC, TimerScope, Threads, CodeVersion, MethodFileMD5, SimulationFileMD5
)]
summary[, MethodOrder := match(Method, method_order)]
setorder(summary, SettingID, MethodOrder)
summary[, MethodOrder := NULL]

# Compute-node hostnames are operational metadata rather than scientific
# provenance, so they are excluded from the public replicate-level aggregate.
if ("Hostname" %in% names(rows)) rows[, Hostname := NULL]

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
raw_file <- file.path(output_dir, "runtime_benchmark_results_current.csv")
summary_file <- file.path(output_dir, "runtime_benchmark_summary_current.csv")
manifest_file <- file.path(output_dir, "runtime_benchmark_manifest.txt")

fwrite(rows, raw_file)
fwrite(summary, summary_file)
writeLines(c(
  paste0("Status: PASS"),
  paste0("Expected tasks: ", nrow(expected)),
  paste0("Observed tasks: ", nrow(rows)),
  paste0("Settings: ", uniqueN(rows$SettingID)),
  paste0("Replicate-level SGE jobs: ", repeats),
  paste0("Replicates per setting/method: ", repeats),
  paste0("Methods: ", paste(method_order, collapse = ", ")),
  paste0("CodeVersion: ", unique(rows$CodeVersion)),
  paste0("MethodFileMD5: ", unique(rows$MethodFileMD5)),
  paste0("SimulationFileMD5: ", unique(rows$SimulationFileMD5)),
  paste0("TimerScope: ", unique(rows$TimerScope)),
  paste0("InputDesign: ", unique(rows$InputDesign)),
  paste0("Threads: ", paste(unique(rows$Threads), collapse = ", "))
), manifest_file)

cat("PASS: aggregated", nrow(rows), "validated runtime task rows.\n")
cat("Raw results:", normalizePath(raw_file), "\n")
cat("Summary:", normalizePath(summary_file), "\n")
