suppressPackageStartupMessages(library(data.table))

script_file <- sub(
  "^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
)
script_dir <- dirname(normalizePath(script_file))
release_root <- dirname(dirname(script_dir))

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(prefix, default) {
  hit <- grep(paste0("^", prefix), args, value = TRUE)
  if (length(hit)) sub(paste0("^", prefix), "", hit[[1]]) else default
}

source_data_dir <- normalizePath(
  arg_value("--source-data-dir=", file.path(release_root, "Source_Data")),
  mustWork = TRUE
)

sim1_file <- file.path(source_data_dir, "simulation_figures_input_alpha005.csv")
sim2_file <- file.path(source_data_dir, "subspace_figures_input_alpha005.csv")
if (!file.exists(sim1_file) || !file.exists(sim2_file)) {
  stop("Missing Simulation 1 or Simulation 2 Source Data in: ", source_data_dir)
}

sim1 <- fread(sim1_file)
sim2 <- fread(sim2_file)

required_sim1 <- c(
  "ParamID", "Replicate", "Method", "ARI_gmm_pc10", "ARI_kmeans_pc10",
  "ARI_louvain_pc10"
)
required_sim2 <- c("ParamID", "Replicate", "Method", "SimilarityPC10")
stopifnot(
  all(required_sim1 %in% names(sim1)),
  all(required_sim2 %in% names(sim2))
)

# ParamID 12 and ParamID 18 are the same shared mean-shift anchor
# (cell_frac = 0.02 and meanlog = 2), displayed in both parameter sweeps.
# Count that scientific condition once in the pooled correlation analysis.
condition_ids <- setdiff(2L:21L, 18L)
sim1 <- sim1[ParamID %in% condition_ids]
sim2 <- sim2[ParamID %in% condition_ids]

rep_counts_1 <- sim1[, .N, by = .(ParamID, Method)]
rep_counts_2 <- sim2[, .N, by = .(ParamID, Method)]
stopifnot(
  nrow(rep_counts_1) == 190L,
  nrow(rep_counts_2) == 190L,
  all(rep_counts_1$N == 10L),
  all(rep_counts_2$N == 10L)
)

clustering <- sim1[, .(
  GMM = mean(ARI_gmm_pc10),
  `k-means` = mean(ARI_kmeans_pc10),
  Louvain = mean(ARI_louvain_pc10)
), by = .(ParamID, Method)]

stability <- sim2[, .(
  PC10_SubspaceSimilarity = mean(SimilarityPC10)
), by = .(ParamID, Method)]

analysis_data <- merge(
  stability, clustering,
  by = c("ParamID", "Method"), all = FALSE, sort = FALSE
)
stopifnot(
  nrow(analysis_data) == 190L,
  !anyDuplicated(analysis_data, by = c("ParamID", "Method"))
)

correlation_row <- function(algorithm) {
  all_methods <- analysis_data
  without_grid <- analysis_data[Method != "PcaGrid"]
  data.table(
    Clustering_algorithm = algorithm,
    Spearman_rho_all_methods = cor(
      all_methods$PC10_SubspaceSimilarity,
      all_methods[[algorithm]],
      method = "spearman"
    ),
    N_all_methods = nrow(all_methods),
    Spearman_rho_excluding_PcaGrid = cor(
      without_grid$PC10_SubspaceSimilarity,
      without_grid[[algorithm]],
      method = "spearman"
    ),
    N_excluding_PcaGrid = nrow(without_grid),
    Analysis_unit = "Paired clustering/stability condition-method mean",
    Replicates_per_mean = 10L,
    Unique_conditions = length(condition_ids),
    Includes_dropout_mid0_identity = TRUE,
    Shared_mean_shift_anchor_counted_once = TRUE
  )
}

correlations <- rbindlist(lapply(c("GMM", "k-means", "Louvain"), correlation_row))
stopifnot(
  identical(round(correlations$Spearman_rho_all_methods, 2), c(0.22, 0.32, 0.62)),
  identical(round(correlations$Spearman_rho_excluding_PcaGrid, 2), c(0.05, 0.19, 0.58)),
  all(correlations$N_all_methods == 190L),
  all(correlations$N_excluding_PcaGrid == 171L)
)

output_file <- file.path(
  source_data_dir,
  "TableS04_subspace_clustering_correlations_source_data.csv"
)
fwrite(correlations, output_file)
cat("Wrote Supplementary Table S4 Source Data: ", output_file, "\n", sep = "")
