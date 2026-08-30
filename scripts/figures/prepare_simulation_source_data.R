#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit) > 1L) stop("Repeated argument: --", name)
  if (length(hit)) substring(hit, nchar(prefix) + 1L) else default
}

sim1_dir <- arg_value("simulation1-dir")
sim2_dir <- arg_value("simulation2-dir")
rds_root <- arg_value("rds-root")
output_dir <- arg_value("output-dir")
qc_source <- arg_value("qc-source", "Simulation aggregates")
if (any(vapply(list(sim1_dir, sim2_dir, rds_root, output_dir), function(x) {
  is.null(x) || !nzchar(x)
}, logical(1)))) {
  stop("Required: --simulation1-dir, --simulation2-dir, --rds-root, --output-dir.")
}
if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE,
                                               no.. = TRUE))) {
  stop("Use an empty output directory; existing Source Data are not overwritten.")
}

read_aggregate <- function(root, name) {
  path <- file.path(root, name)
  if (!file.exists(path)) stop("Missing aggregate: ", path)
  fread(path)
}
publication_method <- function(x) {
  x <- as.character(x)
  x[x == "Grid"] <- "PcaGrid"
  x[x == "Hubert"] <- "PcaHubert"
  x[x == "Tau"] <- "K's tau"
  x
}
method_group <- function(x) {
  ifelse(x %in% c("Winsor", "Quad", "Ball", "Shell", "LR"), "Proposed", "Reference")
}
check_rows <- function(x, param_ids, label) {
  keys <- c("ParamID", "Replicate", "Method", "Alpha")
  if (!all(keys %in% names(x)) || anyDuplicated(x, by = keys)) {
    stop(label, ": missing or duplicate result keys.")
  }
  expected_methods <- c("PCA", "Grid", "Hubert", "PCP", "Tau",
                        "Winsor", "Quad", "Ball", "Shell", "LR")
  if (!setequal(x$ParamID, param_ids) ||
      !setequal(x$Method, expected_methods) ||
      anyNA(x$Replicate) || !setequal(x$Replicate, 1:10)) {
    stop(label, ": unexpected conditions, methods, or replicates.")
  }
  counts <- x[, .N, by = .(ParamID, Method)]
  if (nrow(x) != length(param_ids) * 100L || any(counts$N != 10L)) {
    stop(label, ": expected ten replicates per condition and method.")
  }
}

sim1_all <- read_aggregate(sim1_dir, "simulation1_background_metrics_long.csv")
sim2_all <- read_aggregate(sim2_dir, "simulation2_background_similarity_long.csv")
sim1_summary <- read_aggregate(sim1_dir, "simulation1_background_metrics_summary.csv")
sim2_summary <- read_aggregate(sim2_dir, "simulation2_background_similarity_summary.csv")
sim1 <- sim1_all[is.na(Alpha) | abs(Alpha - 0.05) < 1e-12]
sim2 <- sim2_all[is.na(Alpha) | abs(Alpha - 0.05) < 1e-12]
check_rows(sim1, 1:21, "Clustering simulation")
check_rows(sim2, 2:21, "Stability simulation")
metric_cols <- grep("^(ARI|NMI|k_)", names(sim1), value = TRUE)
if (!length(metric_cols) || any(!is.finite(unlist(sim1[, ..metric_cols]))) ||
    any(!is.finite(unlist(sim2[, .(SimilarityPC10, SimilarityPC20)])))) {
  stop("Non-finite simulation metrics.")
}
sim1[, Method := publication_method(Method)]
sim2[, Method := publication_method(Method)]

# Preserve the published column order, while making provenance paths portable.
sim1[, SourceFile := "results/simulation1/simulation1_background_metrics_long.csv"]
tail1 <- c("GmmItmax", "EvaluationRuntimeSec", "SourceFile", "CodeVersion", "DoubletDonorDesign")
setcolorder(sim1, c(setdiff(names(sim1), tail1), tail1))

relative_rds_paths <- function(paths) {
  root <- paste0(sub("/+$", "", path.expand(rds_root)), "/")
  paths <- as.character(paths)
  if (anyNA(paths) || any(!nzchar(paths))) stop("Missing RDS provenance paths.")
  absolute <- startsWith(paths, "/")
  if (any(absolute & !startsWith(paths, root))) {
    stop("An absolute RDS path is outside --rds-root; supply its original common root.")
  }
  paths[absolute] <- substring(paths[absolute], nchar(root) + 1L)
  if (any(vapply(strsplit(paths, "/", fixed = TRUE), function(x) {
    ".." %in% x
  }, logical(1)))) stop("RDS paths must not traverse outside their root.")
  paths
}
sim2[, `:=`(
  RDS_Path = relative_rds_paths(RDS_Path),
  CleanPath = relative_rds_paths(CleanPath),
  MethodKey = Method,
  SourceFile = "results/simulation2/simulation2_background_similarity_long.csv",
  CleanServerPath = NA_character_,
  LocalPerturbedPath = NA_character_
)]
if ("CleanParamID" %in% names(sim2)) sim2[, CleanParamID := NULL]
tail2 <- c("MethodKey", "SourceFile", "CodeVersion", "DoubletDonorDesign",
           "CleanSeedReference", "CleanGeneHash", "CleanServerPath", "LocalPerturbedPath")
setcolorder(sim2, c(setdiff(names(sim2), tail2), tail2))

pc_sensitivity <- rbindlist(list(
  sim1[, .(ParamID, Replicate, Method, PC = "PC10", ARI = ARI_gmm_pc10)],
  sim1[, .(ParamID, Replicate, Method, PC = "PC20", ARI = ARI_gmm_pc20)]
))
pc_sensitivity[, Group := method_group(Method)]
setorder(pc_sensitivity, ParamID, Replicate, Method, PC)

proposed_methods <- c("Winsor", "Quad", "Ball", "Shell", "LR")
alpha_levels <- c(0.10, 0.20)
sim1_meta <- unique(sim1_all[, .(
  ParamID, Condition, dropout_mid, doublet_frac, cell_frac, gene_frac, meanlog, sdlog
)], by = "ParamID")
make_sim1_alpha <- function(metric, outcome, panel) {
  baseline <- sim1_summary[
    Method %chin% proposed_methods & abs(Alpha - 0.05) < 1e-12,
    .(ParamID, Method, Baseline = get(metric))
  ]
  current <- sim1_summary[
    Method %chin% proposed_methods & Alpha %in% alpha_levels,
    .(ParamID, Method, Alpha, Value = get(metric))
  ]
  current <- merge(current, baseline, by = c("ParamID", "Method"), all.x = TRUE)
  current <- merge(current, sim1_meta, by = "ParamID", all.x = TRUE)
  current[, `:=`(Dataset = "Simulation 1", Panel = panel, Outcome = outcome,
                 Delta = Value - Baseline, Replicates = 10L)]
  current[]
}
make_sim2_alpha <- function(metric, outcome, panel) {
  baseline <- sim2_summary[
    Method %chin% proposed_methods & abs(Alpha - 0.05) < 1e-12,
    .(ParamID, Method, Baseline = get(metric))
  ]
  current <- sim2_summary[
    Method %chin% proposed_methods & Alpha %in% alpha_levels,
    .(ParamID, Condition, Method, Alpha, Value = get(metric))
  ]
  current <- merge(current, baseline, by = c("ParamID", "Method"), all.x = TRUE)
  current[, `:=`(
    Dataset = "Simulation 2", Panel = panel, Outcome = outcome,
    Delta = Value - Baseline, Replicates = 10L,
    dropout_mid = NA_real_, doublet_frac = NA_real_, cell_frac = NA_real_,
    gene_frac = NA_real_, meanlog = NA_real_, sdlog = NA_real_
  )]
  current[]
}
alpha_source <- rbindlist(list(
  make_sim1_alpha("ARI_gmm_pc10", "GMM ARI", "a"),
  make_sim1_alpha("NMI_gmm_pc10", "GMM NMI", "b"),
  make_sim1_alpha("ARI_kmeans_pc10", "k-means ARI", "c"),
  make_sim1_alpha("ARI_louvain_pc10", "Louvain ARI", "c"),
  make_sim2_alpha("SimilarityPC10_mean", "Subspace PC10", "d"),
  make_sim2_alpha("SimilarityPC20_mean", "Subspace PC20", "d")
), use.names = TRUE, fill = TRUE)
setcolorder(alpha_source, c(
  "Dataset", "Panel", "Outcome", "ParamID", "Condition", "Method", "Alpha",
  "Baseline", "Value", "Delta", "Replicates", "dropout_mid", "doublet_frac",
  "cell_frac", "gene_frac", "meanlog", "sdlog"
))
setorder(alpha_source, Panel, Outcome, Method, Alpha, ParamID)
if (nrow(alpha_source) != 1240L ||
    anyDuplicated(alpha_source, by = c("Dataset", "Outcome", "ParamID", "Method", "Alpha")) ||
    any(!is.finite(unlist(alpha_source[, .(Baseline, Value, Delta)])))) {
  stop("Incomplete or duplicate alpha-sensitivity summaries.")
}

grid1 <- sim1[Method == "PcaGrid"]
grid2 <- sim2[Method == "PcaGrid"]
qc <- data.table(
  Item = c(
    "Implementation", "Projection-pursuit scale", "Preprocessing",
    "Internal centering/scaling", "Score construction",
    "Formal Simulation 1 PcaGrid rows", "Formal Simulation 2 PcaGrid rows",
    "Computed principal components", "Maximum loading orthogonality error",
    "Simulation 1 GMM PC10 ARI mean", "Simulation 1 GMM PC10 NMI mean",
    "Simulation 2 PC10 subspace similarity mean"
  ),
  Result = c(
    "rrcov::PcaGrid", "MAD",
    "log-CPM transformation, zero-variance gene removal, and gene-wise standardization",
    "Disabled after manuscript-level preprocessing (center=0, scale=1)",
    "X times the estimated loading matrix", as.character(nrow(grid1)),
    as.character(nrow(grid2)),
    if (all(grid1$grid_computed_k == 20) && all(grid2$grid_computed_k == 20))
      "20 in all checked formal source-data rows" else "See formal source data",
    format(max(c(grid1$grid_orthogonality_error, grid2$grid_orthogonality_error),
               na.rm = TRUE), scientific = TRUE, digits = 3),
    sprintf("%.4f", mean(grid1$ARI_gmm_pc10)),
    sprintf("%.4f", mean(grid1$NMI_gmm_pc10)),
    sprintf("%.4f", mean(grid2$SimilarityPC10))
  ),
  Source = c(rep("Methods implementation", 5), rep(qc_source, 7))
)

outputs <- list(
  simulation_figures_input_alpha005.csv = sim1,
  subspace_figures_input_alpha005.csv = sim2,
  FigureS06_simulation1_pc10_vs_pc20_gmm_ari_source_data.csv = pc_sensitivity,
  FigureS08_alpha_sensitivity_source_data.csv = alpha_source,
  Table_PcaGrid_QC_source_data.csv = qc
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
for (name in names(outputs)) {
  fwrite(outputs[[name]], file.path(output_dir, name), na = "")
  cat(name, ":", nrow(outputs[[name]]), "rows\n")
}
