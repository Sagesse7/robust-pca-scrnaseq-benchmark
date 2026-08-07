suppressPackageStartupMessages({
  library(data.table)
})

script_file <- sub(
  "^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
)
script_dir <- dirname(normalizePath(script_file))
release_root <- dirname(dirname(script_dir))
args <- commandArgs(trailingOnly = TRUE)
source_arg <- grep("^--source-data-dir=", args, value = TRUE)
source_data_dir <- if (length(source_arg)) {
  normalizePath(sub("^--source-data-dir=", "", source_arg[[1]]), mustWork = TRUE)
} else {
  file.path(release_root, "Source_Data")
}

read_source <- function(filename) {
  path <- file.path(source_data_dir, filename)
  if (!file.exists(path)) stop("Missing Source Data file: ", path)
  fread(path, check.names = FALSE)
}

rank_summary <- function(dt, value_col, condition_cols) {
  x <- copy(dt)
  x <- x[is.finite(get(value_col))]
  means <- x[, .(MeanValue = mean(get(value_col))),
             by = c(condition_cols, "Method")]
  means[, Rank := frank(-MeanValue, ties.method = "average"),
        by = condition_cols]
  means[, .(
    mean_value = mean(MeanValue),
    mean_rank = mean(Rank),
    n_comparisons = .N
  ), by = Method]
}

fmt <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

lookup <- function(x, method, column) {
  hit <- x[Method == method, get(column)]
  if (!length(hit)) return(NA_real_)
  hit[[1]]
}

sim1 <- read_source("simulation_figures_input_alpha005.csv")
sim2 <- read_source("subspace_figures_input_alpha005.csv")
real_clustering <- read_source("Figure06_realdata_gmm_pc10_ari_source_data.csv")
real_subspace <- read_source("Figure07_realdata_subspace_pc10_source_data.csv")
runtime <- read_source("Figure08_runtime_scaling_source_data.csv")
simulation_pcwise <- read_source(
  "FigureS04_simulation2_pcwise_pc10_replicate_source_data.csv"
)
summary_table <- read_source("Figure09_method_summary_source_data.csv")

required_methods <- summary_table$Method
stopifnot(length(required_methods) == 10L, !anyDuplicated(required_methods))

# Primary Simulation 1 clustering evidence. Rankings are calculated after
# averaging the 10 independent replicates within each condition and method.
sim_long <- melt(
  sim1,
  id.vars = c("ParamID", "Replicate", "Method"),
  measure.vars = c(
    "ARI_gmm_pc10", "NMI_gmm_pc10",
    "ARI_kmeans_pc10", "NMI_kmeans_pc10",
    "ARI_louvain_pc10", "NMI_louvain_pc10"
  ),
  variable.name = "Measure",
  value.name = "Value"
)
sim_long[, Algorithm := fifelse(
  grepl("_gmm_", Measure), "GMM",
  fifelse(grepl("_kmeans_", Measure), "k-means", "Louvain")
)]
sim_long[, Metric := fifelse(grepl("^ARI_", Measure), "ARI", "NMI")]

doublet_shift <- rank_summary(
  sim_long[Algorithm == "GMM" & ParamID %between% c(7L, 21L)],
  "Value", c("ParamID", "Metric")
)
dropout <- rank_summary(
  sim_long[Algorithm == "GMM" & ParamID %between% c(2L, 6L)],
  "Value", c("ParamID", "Metric")
)

sim_clustering <- rbindlist(lapply(c("GMM", "k-means", "Louvain"), function(a) {
  z <- rank_summary(
    sim_long[Algorithm == a & ParamID %between% c(2L, 21L)],
    "Value", c("ParamID", "Metric")
  )
  z[, Algorithm := a]
  z
}))

# Simulation 2 subspace evidence. The identity reference (dropout.mid = 0,
# ParamID 3) is intentionally excluded from rankings, as in Figure 5.
sim2_ranked <- sim2[ParamID != 3L]
sim_subspace <- rank_summary(
  sim2_ranked, "SimilarityPC10", "ParamID"
)

# Paired Simulation 2 PC-wise evidence used by the manuscript.
simulation_pcwise <- simulation_pcwise[ParamID != 3L]
pcwise <- rank_summary(simulation_pcwise, "PCwisePC10", "ParamID")
real_subspace_summary <- rank_summary(
  real_subspace, "subspace_pc10_mean_cos", "DatasetName"
)

# Real-data GMM evidence. Figure 6 contains three ARI comparisons; duplicated
# rows (when present in the plotting source) are first collapsed by their
# dataset-method-metric key before ranking.
real_clustering_summary <- rank_summary(
  real_clustering, "ARI", c("DatasetName", "Metric")
)

# Runtime evidence from Figure 8. The n = 2000, d = 2000 anchor appears in
# both plotted panels, so it is counted once in the numerical summary below.
runtime_unique <- unique(runtime[, .(
  Method, Cell, Gene, MedianRuntimeSec
)])
runtime_unique <- runtime_unique[is.finite(MedianRuntimeSec)]
runtime_unique[, RuntimeRank := frank(
  MedianRuntimeSec, ties.method = "average"
), by = .(Cell, Gene)]
runtime_summary <- runtime_unique[, .(
  Runtime_median_sec = median(MedianRuntimeSec),
  Runtime_max_sec = max(MedianRuntimeSec),
  Runtime_mean_rank = mean(RuntimeRank),
  Runtime_n_settings = .N
), by = Method]

dimension_order <- c(
  "Doublets / mean shifts", "Dropout", "PC-wise", "Subspace", "Clustering",
  "Runtime"
)

rating_for <- function(method, dimension) {
  raw <- summary_table[Method == method, get(dimension)]
  if (!length(raw)) stop("Missing Figure 9 rating for ", method, " / ", dimension)
  sub("^\\+/-$", "±", raw[[1]])
}

rationale_for <- function(rating, dimension) {
  if (dimension == "Runtime") {
    return(switch(
      rating,
      "+" = paste(
        "The '+' symbol records relatively favorable feasibility in the tested",
        "single-threaded benchmark; it is not a universal speed claim."
      ),
      "±" = paste(
        "The '±' symbol records intermediate runtime in the tested",
        "single-threaded benchmark."
      ),
      "-" = paste(
        "The '-' symbol records computational constraint as scale increased in",
        "the tested single-threaded implementation."
      ),
      stop("Unknown qualitative rating: ", rating)
    ))
  }
  switch(
    rating,
    "+" = paste(
      "The '+' symbol records a relatively favorable tendency in the reported",
      "summaries; it does not claim universal superiority."
    ),
    "±" = paste(
      "The '±' symbol records mixed or condition-dependent evidence; separate",
      "simulation, algorithm, and real-data summaries are retained where applicable."
    ),
    "-" = paste(
      "The '-' symbol records relatively weak performance in the directly assessed",
      "summary and does not imply failure under every setting."
    ),
    stop("Unknown qualitative rating: ", rating)
  )
}

rows <- vector("list", length(required_methods) * length(dimension_order))
ii <- 0L
for (method in required_methods) {
  for (dimension in dimension_order) {
    ii <- ii + 1L
    rating <- rating_for(method, dimension)
    row <- data.table(
      Method = method,
      Dimension = dimension,
      Rating = rating,
      Figure = NA_character_,
      Metric = NA_character_,
      Simulation_mean_value = NA_real_,
      Simulation_mean_rank = NA_real_,
      Simulation_n_comparisons = NA_integer_,
      Kmeans_mean_rank = NA_real_,
      Louvain_mean_rank = NA_real_,
      Real_data_mean_value = NA_real_,
      Real_data_mean_rank = NA_real_,
      Real_data_n_comparisons = NA_integer_,
      Runtime_median_sec = NA_real_,
      Runtime_max_sec = NA_real_,
      Runtime_mean_rank = NA_real_,
      Runtime_n_settings = NA_integer_,
      Major_observation = NA_character_,
      Why_this_supports_the_rating = rationale_for(rating, dimension),
      Evidence_status = "quantitatively derived; qualitative rating curated",
      Interpretation_rule = paste(
        "+ = relatively favorable; ± = mixed or condition-dependent;",
        "- = relatively weak"
      )
    )

    if (dimension == "Doublets / mean shifts") {
      row[, `:=`(
        Figure = "Figure 3",
        Metric = "GMM PC10 ARI/NMI",
        Simulation_mean_value = lookup(doublet_shift, method, "mean_value"),
        Simulation_mean_rank = lookup(doublet_shift, method, "mean_rank"),
        Simulation_n_comparisons = as.integer(lookup(doublet_shift, method, "n_comparisons"))
      )]
      row[, Major_observation := sprintf(
        "GMM mean rank %s across %d condition-metric comparisons (mean score %s).",
        fmt(Simulation_mean_rank), Simulation_n_comparisons,
        fmt(Simulation_mean_value, 3)
      )]
    } else if (dimension == "Dropout") {
      row[, `:=`(
        Figure = "Figure 3; Supplementary Table S3",
        Metric = "GMM PC10 ARI/NMI",
        Simulation_mean_value = lookup(dropout, method, "mean_value"),
        Simulation_mean_rank = lookup(dropout, method, "mean_rank"),
        Simulation_n_comparisons = as.integer(lookup(dropout, method, "n_comparisons"))
      )]
      row[, Major_observation := sprintf(
        "GMM mean rank %s across %d dropout-metric comparisons (mean score %s).",
        fmt(Simulation_mean_rank), Simulation_n_comparisons,
        fmt(Simulation_mean_value, 3)
      )]
    } else if (dimension == "PC-wise") {
      row[, `:=`(
        Figure = "Supplementary Figures S4 and S5",
        Metric = "Mean PC-wise absolute cosine similarity (PC1-PC10)",
        Simulation_mean_value = lookup(pcwise, method, "mean_value"),
        Simulation_mean_rank = lookup(pcwise, method, "mean_rank"),
        Simulation_n_comparisons = as.integer(lookup(pcwise, method, "n_comparisons"))
      )]
      row[, Major_observation := sprintf(
        "Mean score %s and mean rank %s across %d non-identity Simulation 2 conditions.",
        fmt(Simulation_mean_value, 3), fmt(Simulation_mean_rank),
        Simulation_n_comparisons
      )]
    } else if (dimension == "Subspace") {
      row[, `:=`(
        Figure = "Figures 4, 5 and 7",
        Metric = "PC10 mean cosine similarity",
        Simulation_mean_value = lookup(sim_subspace, method, "mean_value"),
        Simulation_mean_rank = lookup(sim_subspace, method, "mean_rank"),
        Simulation_n_comparisons = as.integer(lookup(sim_subspace, method, "n_comparisons")),
        Real_data_mean_value = lookup(real_subspace_summary, method, "mean_value"),
        Real_data_mean_rank = lookup(real_subspace_summary, method, "mean_rank"),
        Real_data_n_comparisons = as.integer(lookup(real_subspace_summary, method, "n_comparisons"))
      )]
      row[, Major_observation := sprintf(
        paste0(
          "Simulation mean rank %s across %d non-identity conditions (mean %s); ",
          "real-data mean rank %s across %d datasets (mean %s)."
        ),
        fmt(Simulation_mean_rank), Simulation_n_comparisons,
        fmt(Simulation_mean_value, 3), fmt(Real_data_mean_rank),
        Real_data_n_comparisons, fmt(Real_data_mean_value, 3)
      )]
    } else if (dimension == "Clustering") {
      gmm <- sim_clustering[Algorithm == "GMM"]
      km <- sim_clustering[Algorithm == "k-means"]
      lv <- sim_clustering[Algorithm == "Louvain"]
      row[, `:=`(
        Figure = "Figures 3 and 6; Supplementary Figures S1-S3",
        Metric = "PC10 ARI/NMI by clustering algorithm; real-data GMM ARI",
        Simulation_mean_value = lookup(gmm, method, "mean_value"),
        Simulation_mean_rank = lookup(gmm, method, "mean_rank"),
        Simulation_n_comparisons = as.integer(lookup(gmm, method, "n_comparisons")),
        Kmeans_mean_rank = lookup(km, method, "mean_rank"),
        Louvain_mean_rank = lookup(lv, method, "mean_rank"),
        Real_data_mean_value = lookup(real_clustering_summary, method, "mean_value"),
        Real_data_mean_rank = lookup(real_clustering_summary, method, "mean_rank"),
        Real_data_n_comparisons = as.integer(lookup(real_clustering_summary, method, "n_comparisons"))
      )]
      row[, Major_observation := sprintf(
        paste0(
          "Simulation mean ranks: GMM %s, k-means %s, Louvain %s ",
          "(%d GMM condition-metric comparisons); real-data GMM mean rank %s ",
          "across %d dataset-metric comparisons."
        ),
        fmt(Simulation_mean_rank), fmt(Kmeans_mean_rank), fmt(Louvain_mean_rank),
        Simulation_n_comparisons, fmt(Real_data_mean_rank),
        Real_data_n_comparisons
      )]
    } else if (dimension == "Runtime") {
      row[, `:=`(
        Figure = "Figure 8",
        Metric = paste(
          "Single-threaded fit plus 20-PC score-construction runtime;",
          "median seconds across five matched replicates"
        ),
        Runtime_median_sec = lookup(runtime_summary, method, "Runtime_median_sec"),
        Runtime_max_sec = lookup(runtime_summary, method, "Runtime_max_sec"),
        Runtime_mean_rank = lookup(runtime_summary, method, "Runtime_mean_rank"),
        Runtime_n_settings = as.integer(lookup(
          runtime_summary, method, "Runtime_n_settings"
        ))
      )]
      row[, Major_observation := sprintf(
        paste0(
          "Median runtime %s s across %d unique scale settings; maximum %s s; ",
          "mean within-setting runtime rank %s (lower is faster)."
        ),
        fmt(Runtime_median_sec, 3), Runtime_n_settings,
        fmt(Runtime_max_sec, 3), fmt(Runtime_mean_rank)
      )]
    }
    rows[[ii]] <- row
  }
}

evidence <- rbindlist(rows, use.names = TRUE, fill = TRUE)
evidence[, Dimension := factor(Dimension, levels = dimension_order)]
evidence[, Method := factor(Method, levels = required_methods)]
setorder(evidence, Dimension, Method)
evidence[, `:=`(Dimension = as.character(Dimension), Method = as.character(Method))]

stopifnot(
  nrow(evidence) == 60L,
  uniqueN(evidence[, .(Method, Dimension)]) == 60L,
  !anyNA(evidence$Rating),
  !anyNA(evidence$Major_observation)
)

output_file <- file.path(source_data_dir, "Figure09_evidence_table.csv")
fwrite(evidence, output_file, na = "")
cat("Figure 9 evidence table written to:", output_file, "\n")
cat("Rows:", nrow(evidence), "(10 methods x 6 dimensions)\n")
