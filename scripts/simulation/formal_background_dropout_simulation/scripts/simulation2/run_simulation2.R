# Simulation 2: clean-versus-perturbed principal-subspace stability.
#
# Clean references are computed once per BaseID and replicate. Perturbed array
# tasks regenerate the paired clean counts only to construct their perturbation;
# they do not repeat the expensive clean dimensionality-reduction fits.

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(readxl)
  library(splatter)
  library(SummarizedExperiment)
})

script_dir <- dirname(normalizePath(sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)[1])))
release_dir <- dirname(dirname(script_dir))
source(file.path(release_dir, "R", "dimensionality_reduction_methods.R"))
source(file.path(release_dir, "R", "simulation_utils.R"))
set_single_thread_blas()

cli <- parse_cli()
param_file <- if (is.null(cli[["param-file"]])) {
  file.path(release_dir, "config", "sim_params_1_v3.xlsx")
} else cli[["param-file"]]
task_index <- if (!is.null(cli[["task-index"]])) {
  as.integer(cli[["task-index"]])
} else if (nzchar(Sys.getenv("SGE_TASK_ID"))) {
  as.integer(Sys.getenv("SGE_TASK_ID"))
} else 1L
mode <- if (is.null(cli$mode)) "perturbed" else {
  match.arg(cli$mode, c("clean-reference", "perturbed"))
}
repeats <- if (is.null(cli$repeats)) 10L else as.integer(cli$repeats)
methods <- cli_csv(cli$methods, all_method_names())
alphas <- cli_numeric_csv(cli$alphas, c(0.05, 0.10, 0.20))
base_seed <- if (is.null(cli[["base-seed"]])) 12345L else as.integer(cli[["base-seed"]])
kpc <- if (is.null(cli$kpc)) 20L else as.integer(cli$kpc)
de_fac_loc <- if (is.null(cli[["de-fac-loc"]])) 0.1 else as.numeric(cli[["de-fac-loc"]])
background_dropout_mid <- if (is.null(cli[["background-dropout-mid"]])) {
  0
} else as.numeric(cli[["background-dropout-mid"]])

run_id <- paste0(
  ifelse(nzchar(Sys.getenv("JOB_ID")), Sys.getenv("JOB_ID"), "manual"),
  ifelse(nzchar(Sys.getenv("SGE_TASK_ID")), paste0("_sge", Sys.getenv("SGE_TASK_ID")), ""),
  "_task", sprintf("%02d", task_index), "_", format(Sys.time(), "%Y%m%d_%H%M%S")
)
default_result_dir <- if (mode == "clean-reference") {
  "simulation2_clean_reference"
} else {
  "simulation2_perturbed"
}
out_dir <- if (is.null(cli[["output-dir"]])) {
  file.path(release_dir, "results", default_result_dir, run_id)
} else cli[["output-dir"]]
nonempty_output_guard(out_dir)

params <- as.data.table(read_excel(param_file, na = "/"))
if (task_index < 1L || task_index > nrow(params)) stop("Invalid task index.")
params[, BaseID := .GRP, by = .(CellType, Cell, Gene, de_prob)]
p <- params[task_index]
requested_condition <- isolated_condition_name(p)
if (mode == "clean-reference" && requested_condition != "clean") {
  stop("Clean-reference mode requires a clean parameter setting.")
}
if (mode == "perturbed" && requested_condition == "clean") {
  stop("Perturbed mode requires a non-clean parameter setting.")
}
run_grid <- as.data.table(method_run_grid(methods, alphas))

inflate_loadings <- function(loadings, keep_names, universe_names, k) {
  out <- matrix(0, nrow = length(universe_names), ncol = k)
  idx <- match(keep_names, universe_names)
  if (anyNA(idx)) stop("Retained analysis genes are absent from the gene universe.")
  out[idx, seq_len(ncol(loadings))] <- loadings
  out
}

manifest <- list()
failures <- list()
mi <- fi <- 0L

checkpoint <- function() {
  if (length(manifest)) {
    fwrite(rbindlist(manifest, fill = TRUE), file.path(out_dir, "simulation2_manifest.csv"))
  }
  if (length(failures)) {
    fwrite(rbindlist(failures, fill = TRUE), file.path(out_dir, "simulation2_runtime_failure.csv"))
  }
}

fwrite(run_grid, file.path(out_dir, "simulation2_method_grid.csv"))
fwrite(params, file.path(out_dir, "simulation_parameter_grid.csv"))

for (replicate_id in seq_len(repeats)) {
  cat(sprintf(
    "Simulation 2 %s task %d/%d, repeat %d/%d\n",
    mode, task_index, nrow(params), replicate_id, repeats
  ))
  clean_seed <- base_seed + p$BaseID * 1000L + replicate_id
  sce_base <- simulate_clean_splatter_base(p, clean_seed, de_fac_loc)
  condition_obj <- make_background_dropout_condition(
    sce_base, p, clean_seed, background_mid = background_dropout_mid
  )
  base_counts <- condition_obj$original_base_counts
  labels <- condition_obj$labels
  condition <- condition_obj$condition
  C <- condition_obj$counts
  contaminated <- condition_obj$contaminated
  condition_seed <- condition_obj$perturbation_seed

  universe <- rownames(C)
  if (is.null(universe)) universe <- paste0("g", seq_len(nrow(C)))
  X <- normalize_log_scale_counts(C)
  if (is.null(colnames(X))) colnames(X) <- universe[seq_len(ncol(X))]
  data_diagnostics <- simulation_condition_diagnostics(
    base_counts, C, condition, contaminated, ncol(X), labels
  )
  data_diagnostics$BackgroundDropoutMid <- condition_obj$background_mid
  data_diagnostics$BackgroundDropoutSeed <- condition_obj$background_seed
  data_diagnostics$BackgroundZeroFraction <- mean(condition_obj$background_counts == 0)
  data_diagnostics$DropoutSweepSeed <- condition_obj$dropout_seed
  data_diagnostics$DoubletDonorDesign <- condition_obj$doublet_donor_design
  data_diagnostics$NoiseDesignDiagnostic <- "background_dropout_mid0"
  geometry <- if (any(run_grid$Method %in% c("Tau", proposed_method_names()))) {
    pairwise_geometry(X)
  } else NULL

  for (j in seq_len(nrow(run_grid))) {
    method <- run_grid$Method[j]
    alpha <- run_grid$Alpha[j]
    tag <- paste0(method, ifelse(is.na(alpha), "", sprintf("_alpha%03d", round(100 * alpha))))
    t0 <- proc.time()[["elapsed"]]
    result <- safe_fit_dimensionality_reduction(
      X, method = method, k = kpc, alpha = ifelse(is.na(alpha), 0.05, alpha),
      geometry = geometry
    )
    elapsed <- proc.time()[["elapsed"]] - t0
    base <- c(as.list(p), list(
      TaskIndex = task_index, Replicate = replicate_id, Condition = condition,
      MasterSeed = base_seed, CleanSeed = clean_seed,
      PerturbationSeed = condition_seed, RunMode = mode,
      CodeVersion = "20260620_doublet_unmodified_singlet_donors",
      NoiseDesign = "background_dropout_mid0", DeFacLoc = de_fac_loc,
      Method = method,
      Calibration = ifelse(method %in% proposed_method_names(), "pairwise_empirical_quantile", "none"),
      Alpha = alpha, RuntimeSec = elapsed
    ), data_diagnostics)
    if (!result$success) {
      fi <- fi + 1L
      failures[[fi]] <- as.data.table(c(
        base, list(success = FALSE, failure_message = result$failure_message)
      ))
      checkpoint()
      next
    }
    if (!is.null(result$fit$cutoff)) {
      base$Q1 <- result$fit$cutoff$Q1
      base$Q2 <- result$fit$cutoff$Q2
      base$Q3 <- result$fit$cutoff$Q3
      base$Q3star <- result$fit$cutoff$Q3star
    }
    base$EigenBackend <- result$fit$backend
    base <- c(base, fit_diagnostic_fields(result$fit))

    loadings <- inflate_loadings(result$fit$loadings, colnames(X), universe, kpc)
    rds_path <- file.path(out_dir, sprintf(
      "BaseID%02d_ParamID%02d_rep%02d_%s_%s.rds",
      p$BaseID, p$ParamID, replicate_id, condition, tag
    ))
    saveRDS(list(
      BaseID = p$BaseID, ParamID = p$ParamID, Replicate = replicate_id,
      Condition = condition, MasterSeed = base_seed, CleanSeed = clean_seed,
      PerturbationSeed = condition_seed, RunMode = mode,
      Method = method, Calibration = result$fit$calibration, Alpha = alpha,
      EigenBackend = result$fit$backend, cutoff = result$fit$cutoff,
      gene_universe = universe, fit_diagnostics = fit_diagnostic_fields(result$fit),
      data_diagnostics = data_diagnostics, loadings_full = loadings,
      eigvals = result$fit$values
    ), rds_path)
    mi <- mi + 1L
    manifest[[mi]] <- as.data.table(c(base, list(
      UniverseP = length(universe), KeptP = ncol(X), GeneHash = digest(universe),
      RDS_Path = rds_path
    )))
    fi <- fi + 1L
    failures[[fi]] <- as.data.table(c(
      base, list(success = TRUE, failure_message = "")
    ))
    checkpoint()
  }
}

checkpoint()
