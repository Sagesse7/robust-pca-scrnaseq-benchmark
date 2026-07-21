suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- args[startsWith(args, prefix)]
  if (length(hit)) sub(prefix, "", hit[[1]], fixed = TRUE) else default
}

input_dirs <- c(
  PBMCs = arg_value("pbmc-dir", ""),
  Pancreas = arg_value("pancreas-dir", ""),
  Bhattacherjee = arg_value("bhattacherjee-dir", "")
)
output_file <- arg_value("output")

if (any(is.na(input_dirs)) || any(!nzchar(input_dirs)) || is.null(output_file)) {
  stop(
    "Required arguments: --pbmc-dir, --pancreas-dir, ",
    "--bhattacherjee-dir, and --output."
  )
}
if (any(!dir.exists(input_dirs))) {
  stop("Missing input directories: ", paste(input_dirs[!dir.exists(input_dirs)], collapse = ", "))
}

publication_method <- function(x) {
  x <- as.character(x)
  x[x == "Grid"] <- "PcaGrid"
  x[x == "Hubert"] <- "PcaHubert"
  x[x == "Tau"] <- "K's tau"
  x
}

method_group <- function(x) {
  ifelse(x %in% c("Winsor", "Quad", "Ball", "Shell", "LR"),
         "Proposed", "Reference")
}

result_status <- c(
  PBMCs = arg_value("pbmc-status", "manuscript analysis"),
  Pancreas = arg_value("pancreas-status", "manuscript analysis"),
  Bhattacherjee = arg_value(
    "bhattacherjee-status",
    "manuscript analysis"
  )
)

required_columns <- list(
  full = c("Replicate", "ClusterRep", "Method", "ARI_gmm_pc10"),
  sample = c(
    "Replicate", "ClusterRep", "Method",
    "gmm_pc10_ARI_full_vs_subset", "ARI_gmm_pc10"
  )
)

derive_dataset <- function(result_dir, dataset_name) {
  full_file <- file.path(result_dir, "real_full_results_repeated_clustering.csv")
  sample_file <- file.path(result_dir, "real_sampling_results_repeated_clustering.csv")
  if (!file.exists(full_file) || !file.exists(sample_file)) {
    stop("Missing repeated-clustering result files in ", result_dir)
  }

  full <- fread(full_file)
  sample <- fread(sample_file)
  if (!all(required_columns$full %in% names(full))) {
    stop("FULL result is missing required columns for ", dataset_name)
  }
  if (!all(required_columns$sample %in% names(sample))) {
    stop("Subset result is missing required columns for ", dataset_name)
  }

  full[, Method := publication_method(Method)]
  sample[, Method := publication_method(Method)]

  out <- rbindlist(list(
    sample[, .(
      DatasetName = dataset_name,
      Method,
      Metric = "FULL vs subset",
      Replicate = as.integer(Replicate),
      ClusterRep = as.integer(ClusterRep),
      ARI = gmm_pc10_ARI_full_vs_subset
    )],
    full[, .(
      DatasetName = dataset_name,
      Method,
      Metric = "FULL vs labels",
      Replicate = NA_integer_,
      ClusterRep = as.integer(ClusterRep),
      ARI = ARI_gmm_pc10
    )],
    sample[, .(
      DatasetName = dataset_name,
      Method,
      Metric = "Subset vs labels",
      Replicate = as.integer(Replicate),
      ClusterRep = as.integer(ClusterRep),
      ARI = ARI_gmm_pc10
    )]
  ), use.names = TRUE)

  out[, `:=`(
    Group = method_group(Method),
    ObservationLevel = "clustering run",
    Figure6SummaryLevel = fifelse(
      Metric == "FULL vs labels",
      "10 clustering runs on one FULL embedding",
      "20 subset means; each mean averages 10 clustering runs"
    ),
    ResultStatus = unname(result_status[[dataset_name]])
  )]
  out
}

figure06 <- rbindlist(Map(derive_dataset, input_dirs, names(input_dirs)), use.names = TRUE)

expected_methods <- c(
  "PCA", "PcaGrid", "PcaHubert", "PCP", "K's tau",
  "Winsor", "Quad", "Ball", "Shell", "LR"
)
if (!setequal(unique(figure06$Method), expected_methods)) {
  stop("Unexpected Figure 6 method set: ", paste(sort(unique(figure06$Method)), collapse = ", "))
}
if (any(!is.finite(figure06$ARI))) stop("Figure 6 contains non-finite ARI values.")

subset_rows <- figure06[Metric != "FULL vs labels"]
subset_counts <- subset_rows[, .N, by = .(DatasetName, Method, Metric, Replicate)]
if (nrow(subset_counts[N != 10L])) {
  stop("Every subset-level value must contain exactly 10 clustering runs.")
}
subset_replicates <- unique(subset_rows[, .(DatasetName, Method, Metric, Replicate)])
replicate_counts <- subset_replicates[, .N, by = .(DatasetName, Method, Metric)]
if (nrow(replicate_counts[N != 20L])) {
  stop("Every subset-dependent method/metric must contain exactly 20 subsets.")
}

full_counts <- figure06[Metric == "FULL vs labels", .N, by = .(DatasetName, Method)]
if (nrow(full_counts[N != 10L])) {
  stop("Every FULL-vs-labels method must contain exactly 10 clustering runs.")
}

key_columns <- c("DatasetName", "Method", "Metric", "Replicate", "ClusterRep")
if (anyDuplicated(figure06, by = key_columns)) {
  stop("Duplicate Figure 6 hierarchical keys detected.")
}

setcolorder(figure06, c(
  "DatasetName", "Method", "Metric", "Replicate", "ClusterRep", "ARI", "Group",
  "ObservationLevel", "Figure6SummaryLevel", "ResultStatus"
))
setorder(figure06, DatasetName, Metric, Method, Replicate, ClusterRep, na.last = TRUE)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
fwrite(figure06, output_file)

summary <- figure06[, .(
  RawRows = .N,
  Subsets = uniqueN(Replicate[!is.na(Replicate)]),
  ClusterRuns = uniqueN(ClusterRep)
), by = .(DatasetName, Metric)]
print(summary)
cat("Wrote hierarchical Figure 6 Source Data to:\n", normalizePath(output_file), "\n")
