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
dir.create(source_data_dir, recursive = TRUE, showWarnings = FALSE)

dt <- read.csv(source_file, stringsAsFactors = FALSE, check.names = FALSE)
num_cols <- c(
  "ParamID", "Replicate", "dropout_mid", "doublet_frac", "cell_frac", "meanlog",
  "ARI_gmm_pc10", "NMI_gmm_pc10"
)
for (cc in intersect(num_cols, names(dt))) dt[[cc]] <- suppressWarnings(as.numeric(dt[[cc]]))

panel_defs <- data.frame(
  Panel = c("dropout", "doublet", "shift_fraction", "shift_factor"),
  ParamID_min = c(2, 7, 12, 17),
  ParamID_max = c(6, 11, 16, 21),
  LevelColumn = c("dropout_mid", "doublet_frac", "cell_frac", "meanlog"),
  PanelLabel = c("Dropout", "Doublet", "Shift fraction", "Shift factor"),
  stringsAsFactors = FALSE
)

panel_data <- lapply(seq_len(nrow(panel_defs)), function(i) {
  p <- panel_defs[i, ]
  x <- dt[dt$ParamID >= p$ParamID_min & dt$ParamID <= p$ParamID_max, , drop = FALSE]
  x$Panel <- p$Panel
  x$PanelLabel <- p$PanelLabel
  x$LevelValue <- x[[p$LevelColumn]]
  x
}) |> bind_rows()

tcrit <- qt(0.975, df = 9)

metric_long <- panel_data |>
  select(Panel, PanelLabel, LevelValue, ParamID, Replicate, Method, ARI_gmm_pc10, NMI_gmm_pc10) |>
  pivot_longer(
    cols = c(ARI_gmm_pc10, NMI_gmm_pc10),
    names_to = "MetricColumn",
    values_to = "Value"
  ) |>
  mutate(Metric = ifelse(MetricColumn == "ARI_gmm_pc10", "ARI", "NMI"))

uncertainty_summary <- metric_long |>
  group_by(Metric, Panel, PanelLabel, LevelValue, ParamID, Method) |>
  summarise(
    n_reps = n(),
    mean = mean(Value, na.rm = TRUE),
    sd = sd(Value, na.rm = TRUE),
    se = sd / sqrt(n_reps),
    ci95_t_low = mean - tcrit * se,
    ci95_t_high = mean + tcrit * se,
    ci95_t_half_width = tcrit * se,
    .groups = "drop"
  )

rank_replicates <- metric_long |>
  group_by(Metric, Panel, PanelLabel, LevelValue, Replicate) |>
  mutate(
    Rank = rank(-Value, ties.method = "average"),
    Top1 = as.integer(Rank <= 1),
    Top3 = as.integer(Rank <= 3)
  ) |>
  ungroup()

rank_stability <- rank_replicates |>
  group_by(Metric, Panel, PanelLabel, LevelValue, ParamID, Method) |>
  summarise(
    n_reps = n(),
    mean_rank = mean(Rank),
    median_rank = median(Rank),
    rank_iqr = IQR(Rank),
    min_rank = min(Rank),
    max_rank = max(Rank),
    top1_freq = mean(Top1),
    top3_freq = mean(Top3),
    .groups = "drop"
  )

panel_uncertainty <- uncertainty_summary |>
  group_by(Metric, Panel, PanelLabel) |>
  summarise(
    median_se = median(se, na.rm = TRUE),
    median_t_ci95_half_width = median(ci95_t_half_width, na.rm = TRUE),
    max_t_ci95_half_width = max(ci95_t_half_width, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    PanelLabel = factor(PanelLabel, levels = panel_defs$PanelLabel)
  ) |>
  arrange(Metric, PanelLabel)

dropout_rank_detail <- uncertainty_summary |>
  filter(Panel == "dropout", Method %in% c("PCP", "K's tau")) |>
  left_join(
    rank_stability,
    by = c("Metric", "Panel", "PanelLabel", "LevelValue", "ParamID", "Method", "n_reps")
  ) |>
  arrange(Metric, LevelValue, Method)

top_dropout <- uncertainty_summary |>
  filter(Panel == "dropout") |>
  group_by(Metric, LevelValue) |>
  slice_max(order_by = mean, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(Metric, LevelValue, top_by_mean_method = Method, top_mean = mean)

dropout_rank_detail <- dropout_rank_detail |>
  left_join(top_dropout, by = c("Metric", "LevelValue"))

paired_dropout <- metric_long |>
  filter(Panel == "dropout", Method %in% c("PCP", "K's tau")) |>
  select(Metric, Panel, PanelLabel, LevelValue, ParamID, Replicate, Method, Value) |>
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
    ci95_t_half_width_diff = tcrit * se_diff,
    ci95_t_low_diff = mean_diff_ktau_minus_pcp - ci95_t_half_width_diff,
    ci95_t_high_diff = mean_diff_ktau_minus_pcp + ci95_t_half_width_diff,
    ktau_gt_pcp_freq = mean(ktau_gt_pcp, na.rm = TRUE),
    pcp_gt_ktau_freq = mean(pcp_gt_ktau, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(Metric, LevelValue)

write.csv(
  uncertainty_summary,
  file.path(source_data_dir, "Appendix_simulation1_gmm_pc10_uncertainty_summary.csv"),
  row.names = FALSE
)
write.csv(
  rank_stability,
  file.path(source_data_dir, "Appendix_simulation1_gmm_pc10_rank_stability.csv"),
  row.names = FALSE
)
write.csv(
  rank_replicates,
  file.path(source_data_dir, "Appendix_simulation1_gmm_pc10_replicate_ranks.csv"),
  row.names = FALSE
)
write.csv(
  panel_uncertainty,
  file.path(source_data_dir, "Appendix_simulation1_gmm_pc10_uncertainty_panel_summary.csv"),
  row.names = FALSE
)
write.csv(
  dropout_rank_detail,
  file.path(source_data_dir, "Appendix_simulation1_dropout_pcp_ktau_rank_detail.csv"),
  row.names = FALSE
)
write.csv(
  paired_dropout_summary,
  file.path(source_data_dir, "Appendix_simulation1_dropout_pcp_ktau_paired_difference.csv"),
  row.names = FALSE
)

message("Wrote Simulation 1 uncertainty and rank-stability Source Data to: ", source_data_dir)
