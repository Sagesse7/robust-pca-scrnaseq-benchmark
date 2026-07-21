# Simulation 1: clustering performance under dropout, doublet, and shift noise.

suppressPackageStartupMessages({
  library(aricode)
  library(data.table)
  library(igraph)
  library(mclust)
  library(readxl)
  library(scran)
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
repeats <- if (is.null(cli$repeats)) 10L else as.integer(cli$repeats)
repeat_ids <- as.integer(cli_numeric_csv(cli[["repeat-ids"]], seq_len(repeats)))
repeat_ids <- unique(repeat_ids)
if (!length(repeat_ids) || any(!is.finite(repeat_ids)) || any(repeat_ids < 1L)) {
  stop("--repeat-ids must contain positive integers.")
}
methods <- cli_csv(cli$methods, all_method_names())
alphas <- cli_numeric_csv(cli$alphas, c(0.05, 0.10, 0.20))
base_seed <- if (is.null(cli[["base-seed"]])) 12345L else as.integer(cli[["base-seed"]])
kpc <- if (is.null(cli$kpc)) 20L else as.integer(cli$kpc)
de_fac_loc <- if (is.null(cli[["de-fac-loc"]])) 0.1 else as.numeric(cli[["de-fac-loc"]])
background_dropout_mid <- if (is.null(cli[["background-dropout-mid"]])) {
  0
} else as.numeric(cli[["background-dropout-mid"]])
clustering_seed_mode <- if (is.null(cli[["clustering-seed-mode"]])) {
  "common"
} else match.arg(cli[["clustering-seed-mode"]], c("legacy", "common"))
gmm_itmax <- if (is.null(cli[["gmm-itmax"]])) 1000L else as.integer(cli[["gmm-itmax"]])
if (!is.finite(gmm_itmax) || gmm_itmax < 1L) stop("--gmm-itmax must be a positive integer.")

run_id <- paste0(
  ifelse(nzchar(Sys.getenv("JOB_ID")), Sys.getenv("JOB_ID"), "manual"),
  ifelse(nzchar(Sys.getenv("SGE_TASK_ID")), paste0("_sge", Sys.getenv("SGE_TASK_ID")), ""),
  "_task", sprintf("%02d", task_index), "_", format(Sys.time(), "%Y%m%d_%H%M%S")
)
out_dir <- if (is.null(cli[["output-dir"]])) {
  file.path(release_dir, "results", "simulation1", run_id)
} else cli[["output-dir"]]
nonempty_output_guard(out_dir)

params <- as.data.table(read_excel(param_file, na = "/"))
if (task_index < 1L || task_index > nrow(params)) stop("Invalid task index.")
params[, BaseID := .GRP, by = .(CellType, Cell, Gene, de_prob)]
p <- params[task_index]
run_grid <- as.data.table(method_run_grid(methods, alphas))

simulate_one <- function(p, seed) {
  sce_base <- simulate_clean_splatter_base(p, seed, de_fac_loc)
  condition <- make_background_dropout_condition(
    sce_base, p, seed, background_mid = background_dropout_mid
  )
  exclude_from_eval <- condition$contaminated & condition$condition == "doublet"
  X <- normalize_log_scale_counts(condition$counts)
  diagnostics <- simulation_condition_diagnostics(
    condition$original_base_counts, condition$counts, condition$condition,
    condition$contaminated, ncol(X), condition$labels
  )
  diagnostics$BackgroundDropoutMid <- condition$background_mid
  diagnostics$BackgroundDropoutSeed <- condition$background_seed
  diagnostics$BackgroundZeroFraction <- mean(condition$background_counts == 0)
  diagnostics$DropoutSweepSeed <- condition$dropout_seed
  diagnostics$DoubletDonorDesign <- condition$doublet_donor_design
  diagnostics$NoiseDesignDiagnostic <- "background_dropout_mid0"

  list(
    X = X,
    labels = condition$labels,
    exclude_from_eval = exclude_from_eval,
    condition = condition$condition,
    diagnostics = diagnostics,
    background_seed = condition$background_seed,
    background_mid = condition$background_mid,
    dropout_seed = condition$dropout_seed,
    perturbation_seed = condition$perturbation_seed
  )
}

evaluate_embedding <- function(
    scores, labels, exclude_from_eval, seed, kpc, method, seed_mode, gmm_itmax) {
  clean <- !exclude_from_eval
  labels_eval <- labels[clean]
  k_true <- uniqueN(labels_eval)
  score10 <- scores[, seq_len(min(10L, ncol(scores))), drop = FALSE]
  score20 <- scores[, seq_len(min(kpc, ncol(scores))), drop = FALSE]

  historical_order <- c("Grid", "PCA", "Winsor", "Quad", "Ball", "Shell", "LR", "Tau", "Hubert", "PCP")
  method_offset <- if (seed_mode == "legacy") match(method, historical_order) else 0L
  km10_values <- km20_values <- matrix(NA_real_, nrow = 20L, ncol = 2L)
  for (r in seq_len(20L)) {
    set.seed(seed + method_offset + r)
    cl10 <- kmeans(score10, centers = k_true, nstart = 10L)$cluster[clean]
    cl20 <- kmeans(score20, centers = k_true, nstart = 10L)$cluster[clean]
    km10_values[r, ] <- c(adjustedRandIndex(cl10, labels_eval), aricode::NMI(cl10, labels_eval))
    km20_values[r, ] <- c(adjustedRandIndex(cl20, labels_eval), aricode::NMI(cl20, labels_eval))
  }
  km10 <- colMeans(km10_values)
  km20 <- colMeans(km20_values)

  lv_eval <- function(Z) {
    Z <- scale(Z)
    graph <- scran::buildSNNGraph(t(Z), k = min(15L, nrow(Z) - 1L), type = "jaccard")
    cl <- as.integer(cluster_louvain(graph, weights = E(graph)$weight)$membership)
    c(
      ARI = adjustedRandIndex(cl[clean], labels_eval),
      NMI = aricode::NMI(cl[clean], labels_eval),
      K = uniqueN(cl)
    )
  }
  gmm_eval <- function(Z) {
    fit <- Mclust(
      Z, G = 1:15, verbose = FALSE,
      control = emControl(itmax = c(gmm_itmax, gmm_itmax))
    )
    cl <- as.integer(fit$classification)
    c(
      ARI = adjustedRandIndex(cl[clean], labels_eval),
      NMI = aricode::NMI(cl[clean], labels_eval),
      K = fit$G
    )
  }

  lv10 <- lv_eval(score10)
  lv20 <- lv_eval(score20)
  gm10 <- gmm_eval(score10)
  gm20 <- gmm_eval(score20)
  list(
    ARI_kmeans_pc10 = km10[1], NMI_kmeans_pc10 = km10[2],
    ARI_louvain_pc10 = lv10["ARI"], NMI_louvain_pc10 = lv10["NMI"], k_louvain_pc10 = lv10["K"],
    ARI_gmm_pc10 = gm10["ARI"], NMI_gmm_pc10 = gm10["NMI"], k_gmm_pc10 = gm10["K"],
    ARI_kmeans_pc20 = km20[1], NMI_kmeans_pc20 = km20[2],
    ARI_louvain_pc20 = lv20["ARI"], NMI_louvain_pc20 = lv20["NMI"], k_louvain_pc20 = lv20["K"],
    ARI_gmm_pc20 = gm20["ARI"], NMI_gmm_pc20 = gm20["NMI"], k_gmm_pc20 = gm20["K"]
  )
}

metrics <- list()
failures <- list()
mi <- fi <- 0L

checkpoint <- function() {
  if (length(metrics)) {
    fwrite(rbindlist(metrics, fill = TRUE), file.path(out_dir, "simulation1_metrics.csv"))
  }
  if (length(failures)) {
    fwrite(rbindlist(failures, fill = TRUE), file.path(out_dir, "simulation1_runtime_failure.csv"))
  }
}

fwrite(run_grid, file.path(out_dir, "simulation1_method_grid.csv"))
fwrite(params, file.path(out_dir, "simulation_parameter_grid.csv"))

log_progress <- function(...) {
  cat(sprintf(...), "\n")
  flush.console()
}

for (replicate_id in repeat_ids) {
  log_progress("Simulation 1 task %d/%d, repeat %d", task_index, nrow(params), replicate_id)
  data_seed <- base_seed + p$BaseID * 1000L + replicate_id
  log_progress("[repeat %d] START data generation (seed=%d)", replicate_id, data_seed)
  dat <- simulate_one(p, data_seed)
  log_progress("[repeat %d] END data generation", replicate_id)
  geometry <- if (any(run_grid$Method %in% c("Tau", proposed_method_names()))) {
    log_progress("[repeat %d] START pairwise geometry", replicate_id)
    pairwise_geometry(dat$X)
  } else NULL
  if (!is.null(geometry)) log_progress("[repeat %d] END pairwise geometry", replicate_id)

  for (j in seq_len(nrow(run_grid))) {
    method <- run_grid$Method[j]
    alpha <- run_grid$Alpha[j]
    method_tag <- paste0(method, ifelse(is.na(alpha), "", sprintf(" alpha=%.2f", alpha)))
    log_progress("[repeat %d] START fit: %s", replicate_id, method_tag)
    t0 <- proc.time()[["elapsed"]]
    result <- safe_fit_dimensionality_reduction(
      dat$X, method = method, k = kpc, alpha = ifelse(is.na(alpha), 0.05, alpha),
      geometry = geometry
    )
    elapsed <- proc.time()[["elapsed"]] - t0
    log_progress("[repeat %d] END fit: %s (%.1f sec)", replicate_id, method_tag, elapsed)
    base <- c(as.list(p), list(
      TaskIndex = task_index, Replicate = replicate_id, Method = method,
      MasterSeed = base_seed, CleanSeed = data_seed,
      PerturbationSeed = dat$perturbation_seed,
      CodeVersion = "20260620_doublet_unmodified_singlet_donors",
      NoiseDesign = "background_dropout_mid0", Condition = dat$condition, DeFacLoc = de_fac_loc,
      ClusteringSeedMode = clustering_seed_mode, GmmItmax = gmm_itmax,
      Calibration = ifelse(method %in% proposed_method_names(), "pairwise_empirical_quantile", "none"),
      Alpha = alpha, RuntimeSec = elapsed
    ), dat$diagnostics)
    if (!result$success) {
      fi <- fi + 1L
      failures[[fi]] <- as.data.table(c(base, list(success = FALSE, failure_message = result$failure_message)))
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
    log_progress("[repeat %d] START clustering evaluation: %s", replicate_id, method_tag)
    eval_t0 <- proc.time()[["elapsed"]]
    ev_result <- tryCatch({
      scores <- dat$X %*% result$fit$loadings
      evaluate_embedding(
        scores, dat$labels, dat$exclude_from_eval, data_seed, kpc, method,
        clustering_seed_mode, gmm_itmax
      )
    }, error = function(e) e)
    base$EvaluationRuntimeSec <- proc.time()[["elapsed"]] - eval_t0
    if (inherits(ev_result, "error")) {
      fi <- fi + 1L
      failures[[fi]] <- as.data.table(c(
        base, list(success = FALSE, failure_message = conditionMessage(ev_result))
      ))
      checkpoint()
      next
    }
    log_progress(
      "[repeat %d] END clustering evaluation: %s (%.1f sec)",
      replicate_id, method_tag, base$EvaluationRuntimeSec
    )
    ev <- ev_result
    mi <- mi + 1L
    metrics[[mi]] <- as.data.table(c(base, ev))
    fi <- fi + 1L
    failures[[fi]] <- as.data.table(c(base, list(success = TRUE, failure_message = "")))
    checkpoint()
  }
}

checkpoint()
