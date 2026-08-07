suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

script_file <- sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)[1])
script_dir <- dirname(normalizePath(script_file))
release_root <- dirname(dirname(script_dir))

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1]]) else default
}

source_data_dir <- file.path(release_root, "Source_Data")
source_file <- arg_value(
  "source-file",
  file.path(source_data_dir, "simulation_figures_input_alpha005.csv")
)
output_file <- arg_value(
  "output-file",
  file.path(source_data_dir, "Appendix_simulation1_dropout_pcp_ktau_paired_difference.csv")
)
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

dt <- read.csv(source_file, stringsAsFactors = FALSE, check.names = FALSE)
required_columns <- c(
  "ParamID", "Replicate", "dropout_mid", "Method",
  "ARI_gmm_pc10", "NMI_gmm_pc10"
)
missing_columns <- setdiff(required_columns, names(dt))
if (length(missing_columns)) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

numeric_columns <- setdiff(required_columns, "Method")
for (column in numeric_columns) {
  dt[[column]] <- suppressWarnings(as.numeric(dt[[column]]))
}

paired_dropout <- dt |>
  filter(
    ParamID %in% 2:6,
    Method %in% c("PCP", "K's tau")
  ) |>
  select(
    ParamID, Replicate, dropout_mid, Method,
    ARI_gmm_pc10, NMI_gmm_pc10
  ) |>
  pivot_longer(
    cols = c(ARI_gmm_pc10, NMI_gmm_pc10),
    names_to = "MetricColumn",
    values_to = "Value"
  ) |>
  mutate(
    Metric = ifelse(MetricColumn == "ARI_gmm_pc10", "ARI", "NMI"),
    Panel = "dropout",
    PanelLabel = "Dropout",
    LevelValue = dropout_mid
  ) |>
  select(
    Metric, Panel, PanelLabel, LevelValue,
    ParamID, Replicate, Method, Value
  ) |>
  pivot_wider(names_from = Method, values_from = Value) |>
  mutate(
    ktau_minus_pcp = `K's tau` - PCP,
    ktau_gt_pcp = as.integer(ktau_minus_pcp > 0),
    pcp_gt_ktau = as.integer(ktau_minus_pcp < 0)
  )

paired_dropout_summary <- paired_dropout |>
  group_by(Metric, Panel, PanelLabel, LevelValue, ParamID) |>
  summarise(
    n_reps = n(),
    mean_diff_ktau_minus_pcp = mean(ktau_minus_pcp, na.rm = TRUE),
    sd_diff = sd(ktau_minus_pcp, na.rm = TRUE),
    se_diff = sd_diff / sqrt(n_reps),
    ci95_t_half_width_diff = qt(0.975, df = n_reps - 1) * se_diff,
    ci95_t_low_diff = mean_diff_ktau_minus_pcp - ci95_t_half_width_diff,
    ci95_t_high_diff = mean_diff_ktau_minus_pcp + ci95_t_half_width_diff,
    ktau_gt_pcp_freq = mean(ktau_gt_pcp, na.rm = TRUE),
    pcp_gt_ktau_freq = mean(pcp_gt_ktau, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(Metric, LevelValue)

if (nrow(paired_dropout_summary) != 10L ||
    any(paired_dropout_summary$n_reps != 10L)) {
  stop("Unexpected dropout paired-comparison dimensions.")
}

write.csv(paired_dropout_summary, output_file, row.names = FALSE)
message("Wrote Supplementary Table S3 Source Data to: ", output_file)
