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

as_flag <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  tolower(x) %in% c("1", "true", "yes", "y")
}

package_version_or_na <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(pkg))
}

index_signature <- function(idx) {
  idx <- as.numeric(idx)
  paste(
    length(idx),
    format(sum(idx), scientific = FALSE, trim = TRUE),
    format(sum(idx^2), scientific = FALSE, trim = TRUE),
    sep = ":"
  )
}

elapsed_time <- function() proc.time()[["elapsed"]]

cli <- parse_cli(commandArgs(trailingOnly = TRUE))
grid_file <- if (is.null(cli[["grid-file"]])) {
  file.path(bundle_root, "config", "runtime_grid.csv")
} else cli[["grid-file"]]
output_dir <- if (is.null(cli[["output-dir"]])) {
  file.path(bundle_root, "results", "tasks")
} else cli[["output-dir"]]
repeats <- if (is.null(cli$repeats)) 5L else as.integer(cli$repeats)
base_seed <- if (is.null(cli[["base-seed"]])) 12345L else as.integer(cli[["base-seed"]])
alpha <- if (is.null(cli$alpha)) 0.05 else as.numeric(cli$alpha)
overwrite <- as_flag(cli$overwrite, FALSE)

if (!file.exists(grid_file)) stop("Runtime grid not found: ", grid_file)
if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) stop("alpha must be in (0, 1).")
if (!is.finite(repeats) || repeats < 1L) stop("repeats must be positive.")

method_order <- c(
  "PCA", "Grid", "Hubert", "PCP", "Tau",
  "Winsor", "Quad", "Ball", "Shell", "LR"
)
proposed_methods <- c("Winsor", "Quad", "Ball", "Shell", "LR")
generalized_methods <- c("Tau", proposed_methods)
method_label <- c(
  PCA = "PCA", Grid = "PcaGrid", Hubert = "PcaHubert", PCP = "PCP",
  Tau = "K's tau", Winsor = "Winsor", Quad = "Quad", Ball = "Ball",
  Shell = "Shell", LR = "LR"
)

grid <- fread(grid_file)
required_grid_columns <- c(
  "SettingID", "Panel", "Cell", "Gene", "CellType", "de_prob",
  "DeFacLoc", "BackgroundDropoutMid"
)
missing_grid_columns <- setdiff(required_grid_columns, names(grid))
if (length(missing_grid_columns)) {
  stop("Runtime grid is missing columns: ", paste(missing_grid_columns, collapse = ", "))
}
if (anyDuplicated(grid$SettingID)) stop("SettingID must be unique in the runtime grid.")
if (!identical(sort(unique(grid$Cell)), c(500L, 1000L, 2000L, 3000L)) ||
    !identical(sort(unique(grid$Gene)), c(500L, 1000L, 2000L, 3000L))) {
  stop("The matched runtime design requires Cell and Gene values 500, 1000, 2000, and 3000.")
}

task_map <- as.data.table(expand.grid(
  SettingID = grid$SettingID,
  Replicate = seq_len(repeats),
  Method = method_order,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
))
task_map[, TaskID := .I]

explicit_mode <- all(c("setting-id", "replicate", "method") %in% names(cli))
if (explicit_mode) {
  setting_id <- as.integer(cli[["setting-id"]])
  replicate_id <- as.integer(cli$replicate)
  method <- cli$method
  hit <- task_map[
    SettingID == setting_id & Replicate == replicate_id & Method == method
  ]
  if (nrow(hit) != 1L) stop("Explicit setting/replicate/method did not identify one task.")
  task <- hit
} else {
  task_id <- if (!is.null(cli[["task-id"]])) {
    as.integer(cli[["task-id"]])
  } else if (nzchar(Sys.getenv("SGE_TASK_ID"))) {
    as.integer(Sys.getenv("SGE_TASK_ID"))
  } else {
    stop("Provide --task-id or --setting-id, --replicate, and --method.")
  }
  if (!is.finite(task_id) || task_id < 1L || task_id > nrow(task_map)) {
    stop("task-id must be between 1 and ", nrow(task_map), ".")
  }
  task <- task_map[TaskID == task_id]
}

setting <- grid[SettingID == task$SettingID]
stopifnot(nrow(setting) == 1L)

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1L)
  RhpcBLASctl::omp_set_num_threads(1L)
}
data.table::setDTthreads(1L)

simulation_bundle_root <- normalizePath(
  file.path(
    bundle_root, "..", "..", "simulation",
    "formal_background_dropout_simulation"
  ),
  mustWork = TRUE
)
methods_file <- file.path(simulation_bundle_root, "R", "dimensionality_reduction_methods.R")
simulation_file <- file.path(simulation_bundle_root, "R", "simulation_utils.R")
if (!file.exists(methods_file) || !file.exists(simulation_file)) {
  stop("The shared method/simulation source files are missing.")
}
source(simulation_file)
source(methods_file)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(output_dir, sprintf("runtime_task_%04d.csv", task$TaskID))
if (file.exists(out_file) && !overwrite) {
  stop("Output already exists; use --overwrite=true only after checking provenance: ", out_file)
}

# A single 3000-cell x 3000-gene Splatter draw is reconstructed deterministically
# in every method process. Each target input is an ordered, nested subset of this
# same replicate-level draw. This makes dimensions paired within a replicate while
# retaining independent R processes and no shared fitted object across methods.
master_cell <- 3000L
master_gene <- 3000L
master_seed <- as.integer(base_seed + task$Replicate * 100000L)
cell_subset_seed <- as.integer(master_seed + 1001L)
gene_subset_seed <- as.integer(master_seed + 2001L)
method_seed <- as.integer(master_seed + match(task$Method, method_order) * 1000L)

master_parameters <- as.list(setting)
master_parameters$Cell <- master_cell
master_parameters$Gene <- master_gene
master_parameters$dropout_mid <- NA_real_
master_parameters$doublet_frac <- NA_real_
master_parameters$cell_frac <- NA_real_
master_parameters$gene_frac <- NA_real_
master_parameters$meanlog <- NA_real_
master_parameters$sdlog <- NA_real_

set.seed(master_seed)
sce_master <- simulate_clean_splatter_base(
  master_parameters, master_seed, setting$DeFacLoc
)
master_condition <- make_background_dropout_condition(
  sce_master, master_parameters, master_seed,
  background_mid = setting$BackgroundDropoutMid
)
C_master <- as.matrix(master_condition$background_counts)
if (nrow(C_master) != master_gene || ncol(C_master) != master_cell) {
  stop("The master count matrix has unexpected dimensions.")
}

set.seed(cell_subset_seed)
cell_order <- sample.int(ncol(C_master))
set.seed(gene_subset_seed)
gene_order <- sample.int(nrow(C_master))
selected_cells <- cell_order[seq_len(setting$Cell)]
selected_genes <- gene_order[seq_len(setting$Gene)]
C <- C_master[selected_genes, selected_cells, drop = FALSE]
X <- normalize_log_scale_counts(C)
if (nrow(X) != setting$Cell) stop("Unexpected number of cells after preprocessing.")
if (ncol(X) < 20L) stop("Fewer than 20 genes remained after preprocessing.")

fit_generalized_profiled <- function(X, method, k, alpha) {
  t <- elapsed_time()
  geometry <- pairwise_geometry(X)
  geometry_sec <- elapsed_time() - t

  t <- elapsed_time()
  cutoff <- if (method == "Tau") NULL else pairwise_empirical_cutoffs(geometry, alpha)
  W <- squared_pairwise_weights(method, geometry, cutoff)
  weight_sec <- elapsed_time() - t

  t <- elapsed_time()
  K <- weighted_pairwise_matrix(X, W)
  scatter_sec <- elapsed_time() - t

  t <- elapsed_time()
  fit <- top_symmetric_eigen(K, k)
  eigen_sec <- elapsed_time() - t

  fit$method <- method
  fit$calibration <- if (method == "Tau") "none" else cutoff$calibration
  fit$alpha <- if (method == "Tau") NA_real_ else alpha
  fit$cutoff <- cutoff
  list(
    fit = fit,
    GeometrySec = geometry_sec,
    WeightSec = weight_sec,
    ScatterSec = scatter_sec,
    EigenSec = eigen_sec,
    FitSec = geometry_sec + weight_sec + scatter_sec + eigen_sec
  )
}

invisible(gc(full = TRUE))
set.seed(method_seed)
fit_result <- tryCatch(
  {
    if (task$Method %in% generalized_methods) {
      profiled <- fit_generalized_profiled(X, task$Method, 20L, alpha)
      fit <- profiled$fit
      timings <- profiled[names(profiled) != "fit"]
    } else {
      t_fit <- elapsed_time()
      fit <- fit_dimensionality_reduction(
        X,
        method = task$Method,
        k = 20L,
        alpha = alpha,
        geometry = NULL,
        return_weights = FALSE
      )
      timings <- list(
        GeometrySec = NA_real_, WeightSec = NA_real_, ScatterSec = NA_real_,
        EigenSec = NA_real_, FitSec = elapsed_time() - t_fit
      )
    }

    t_score <- elapsed_time()
    scores <- X %*% fit$loadings
    score_sec <- elapsed_time() - t_score
    if (ncol(scores) < 20L || any(!is.finite(scores))) {
      stop("Invalid score matrix returned.")
    }
    list(
      success = TRUE, fit = fit, scores = scores, failure_message = "",
      timings = c(timings, list(ScoreSec = score_sec))
    )
  },
  error = function(e) {
    list(
      success = FALSE, fit = NULL, scores = NULL,
      failure_message = conditionMessage(e),
      timings = list(
        GeometrySec = NA_real_, WeightSec = NA_real_, ScatterSec = NA_real_,
        EigenSec = NA_real_, FitSec = NA_real_, ScoreSec = NA_real_
      )
    )
  }
)

fit <- fit_result$fit
scores <- fit_result$scores
timings <- fit_result$timings
runtime_sec <- if (fit_result$success) {
  as.numeric(timings$FitSec + timings$ScoreSec)
} else {
  NA_real_
}
orthogonality_error <- if (fit_result$success) {
  max(abs(crossprod(fit$loadings) - diag(ncol(fit$loadings))))
} else NA_real_

cutoff_value <- function(name) {
  if (!fit_result$success || is.null(fit$cutoff) || is.null(fit$cutoff[[name]])) {
    return(NA_real_)
  }
  as.numeric(fit$cutoff[[name]])
}

row <- data.table(
  TaskID = task$TaskID,
  SettingID = setting$SettingID,
  Panel = setting$Panel,
  Replicate = task$Replicate,
  Cell = setting$Cell,
  Gene = setting$Gene,
  AnalysisGenes = ncol(X),
  CellType = setting$CellType,
  de_prob = setting$de_prob,
  DeFacLoc = setting$DeFacLoc,
  BackgroundDropoutMid = setting$BackgroundDropoutMid,
  InputDesign = "nested_subsets_from_common_3000x3000_splatter_background_mid0",
  MasterCell = master_cell,
  MasterGene = master_gene,
  MasterSeed = master_seed,
  CellSubsetSeed = cell_subset_seed,
  GeneSubsetSeed = gene_subset_seed,
  CellSubsetSignature = index_signature(selected_cells),
  GeneSubsetSignature = index_signature(selected_genes),
  Method = task$Method,
  MethodLabel = unname(method_label[[task$Method]]),
  Calibration = if (task$Method %in% proposed_methods) {
    "pairwise_empirical_quantile"
  } else {
    "none"
  },
  Alpha = if (task$Method %in% proposed_methods) alpha else NA_real_,
  PC = 20L,
  TimerScope = "fit_plus_score_excludes_data_generation_preprocessing; generalized_components_profiled",
  RuntimeSec = runtime_sec,
  FitSec = as.numeric(timings$FitSec),
  GeometrySec = as.numeric(timings$GeometrySec),
  WeightSec = as.numeric(timings$WeightSec),
  ScatterSec = as.numeric(timings$ScatterSec),
  EigenSec = as.numeric(timings$EigenSec),
  ScoreSec = as.numeric(timings$ScoreSec),
  Success = fit_result$success,
  FailureMessage = fit_result$failure_message,
  BaseSeed = base_seed,
  MethodSeed = method_seed,
  Threads = 1L,
  EigenBackend = if (fit_result$success) fit$backend else NA_character_,
  ComputedPCs = if (fit_result$success) ncol(fit$loadings) else NA_integer_,
  LoadingOrthogonalityError = orthogonality_error,
  InputSumSquares = sum(X^2),
  ScoreSum = if (fit_result$success) sum(scores) else NA_real_,
  ScoreSumSquares = if (fit_result$success) sum(scores^2) else NA_real_,
  Q1 = cutoff_value("Q1"),
  Q2 = cutoff_value("Q2"),
  Q3 = cutoff_value("Q3"),
  Q3star = cutoff_value("Q3star"),
  CodeVersion = "20260712_matched_nested_runtime_v2",
  MethodFileMD5 = unname(tools::md5sum(methods_file)),
  SimulationFileMD5 = unname(tools::md5sum(simulation_file)),
  RVersion = R.version.string,
  SplatterVersion = package_version_or_na("splatter"),
  RrcovVersion = package_version_or_na("rrcov"),
  RpcaVersion = package_version_or_na("rpca"),
  RSpectraVersion = package_version_or_na("RSpectra"),
  Hostname = Sys.info()[["nodename"]],
  TimestampUTC = format(Sys.time(), tz = "UTC", usetz = TRUE)
)

tmp_file <- paste0(out_file, ".tmp_", Sys.getpid())
fwrite(row, tmp_file)
if (!file.rename(tmp_file, out_file)) {
  unlink(tmp_file)
  stop("Atomic rename failed for: ", out_file)
}

cat(sprintf(
  "Completed runtime task %d/%d: setting=%d replicate=%d method=%s runtime=%.3f s success=%s\n",
  task$TaskID, nrow(task_map), setting$SettingID, task$Replicate,
  task$Method, runtime_sec, fit_result$success
))
