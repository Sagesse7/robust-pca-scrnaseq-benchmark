suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(rsvg)
})

args <- commandArgs(trailingOnly = TRUE)
only_arg <- grep("^--only=", args, value = TRUE)
figure_only <- if (length(only_arg)) sub("^--only=", "", only_arg[[1]]) else "all"
if (!figure_only %in% c("all", "Figure06", "Figure06_07")) {
  stop("--only must be one of all, Figure06, or Figure06_07.")
}

script_file <- sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)[1])
script_dir <- dirname(normalizePath(script_file))
release_root <- dirname(dirname(script_dir))
source(file.path(release_root, "scripts", "utils", "figure_style.R"))

source_data_dir <- file.path(release_root, "Source_Data")
fig_path <- function(folder, stem) file.path(release_root, "figures", folder, stem)
read_source <- function(filename) {
  path <- file.path(source_data_dir, filename)
  if (!file.exists(path)) stop("Missing Source Data: ", path)
  fread(path)
}

datasets <- c("PBMCs", "Pancreas", "Bhattacherjee")
prepare_methods <- function(x) {
  x[, `:=`(
    DatasetName = factor(DatasetName, levels = datasets),
    Method = clean_method(Method)
  )]
  x
}

metric_data <- prepare_methods(read_source("Figure06_realdata_gmm_pc10_ari_source_data.csv"))
required_figure06_columns <- c("Replicate", "ClusterRep", "ARI")
if (!all(required_figure06_columns %in% names(metric_data))) {
  stop(
    "Figure 6 Source Data must retain Replicate and ClusterRep so that ",
    "clustering runs can be summarized within each sampled subset."
  )
}
metric_data[, Metric := factor(
  Metric, levels = c("FULL vs subset", "FULL vs labels", "Subset vs labels")
)]
metric_data[, Group := method_group(Method)]

subset_level <- metric_data[
  is.finite(ARI) & as.character(Metric) != "FULL vs labels",
  .(ARI = mean(ARI), ClusterRuns = .N),
  by = .(DatasetName, Method, Metric, Group, Replicate)
]
if (nrow(subset_level[ClusterRuns != 10L])) {
  stop("Each Figure 6 subset mean must contain exactly 10 clustering runs.")
}
subset_counts <- subset_level[, .N, by = .(DatasetName, Method, Metric)]
if (nrow(subset_counts[N != 20L])) {
  stop("Each Figure 6 subset-dependent method/metric must contain 20 subsets.")
}

full_runs <- metric_data[
  is.finite(ARI) & as.character(Metric) == "FULL vs labels"
]
full_counts <- full_runs[, .N, by = .(DatasetName, Method)]
if (nrow(full_counts[N != 10L])) {
  stop("Each Figure 6 FULL-vs-labels method must contain 10 clustering runs.")
}
full_summary <- full_runs[, .(
  Mean = mean(ARI),
  Minimum = min(ARI),
  Maximum = max(ARI)
), by = .(DatasetName, Method, Metric, Group)]

metric_labels <- c(
  "FULL vs subset" = "FULL vs subset\n20 subset means",
  "FULL vs labels" = "FULL vs labels\n10 clustering runs",
  "Subset vs labels" = "Subset vs labels\n20 subset means"
)

fig06 <- ggplot(mapping = aes(Method, fill = Group)) +
  geom_boxplot(
    data = subset_level, aes(y = ARI),
    outlier.shape = NA, width = 0.64, linewidth = 0.32
  ) +
  geom_point(
    data = subset_level, aes(y = ARI, colour = Group),
    position = position_jitter(width = 0.10, height = 0, seed = 601L),
    alpha = 0.48, size = 0.46
  ) +
  geom_linerange(
    data = full_summary,
    aes(ymin = Minimum, ymax = Maximum, colour = Group),
    linewidth = 0.42
  ) +
  geom_crossbar(
    data = full_summary,
    aes(y = Mean, ymin = Mean, ymax = Mean),
    width = 0.46, linewidth = 0.42, colour = "#202020"
  ) +
  geom_point(
    data = full_runs, aes(y = ARI, colour = Group),
    position = position_jitter(width = 0.14, height = 0, seed = 602L),
    shape = 21, stroke = 0.20, alpha = 0.82, size = 0.72,
    show.legend = FALSE
  ) +
  facet_grid(
    Metric ~ DatasetName,
    labeller = labeller(Metric = as_labeller(metric_labels))
  ) +
  scale_x_discrete(labels = METHOD_LABELS) +
  scale_fill_manual(
    values = METHOD_GROUP_FILLS,
    breaks = c("Reference", "Proposed"),
    labels = c("Reference methods", "Proposed methods")
  ) +
  scale_colour_manual(
    values = METHOD_GROUP_POINTS,
    breaks = c("Reference", "Proposed"),
    labels = c("Reference methods", "Proposed methods")
  ) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  coord_cartesian(ylim = c(-0.05, 1.02)) +
  labs(x = NULL, y = "Adjusted Rand index", fill = NULL) +
  guides(colour = "none") +
  theme_journal(6.2) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.text = element_text(size = 6.2),
    legend.key.width = unit(5.0, "mm"),
    legend.spacing.x = unit(1.2, "mm"),
    panel.spacing = unit(1.4, "mm")
  )
save_pub_r(fig06, fig_path("Real_data", "Figure06_realdata_gmm_pc10_ari"), 183, 135)

if (figure_only == "Figure06") {
  cat("Regenerated Figure 6 from hierarchical Source Data.\n")
  quit(save = "no", status = 0L)
}

stability <- prepare_methods(read_source("Figure07_realdata_subspace_pc10_source_data.csv"))
stability[, Group := method_group(Method)]
fig07 <- ggplot(
  stability[is.finite(subspace_pc10_mean_cos)],
  aes(Method, subspace_pc10_mean_cos, fill = Group)
) +
  geom_boxplot(outlier.shape = NA, width = 0.64, linewidth = 0.32) +
  geom_point(
    aes(colour = Group),
    position = position_jitter(width = 0.10, height = 0, seed = 701L),
    alpha = 0.42, size = 0.42
  ) +
  facet_wrap(~DatasetName, nrow = 1) +
  scale_x_discrete(labels = METHOD_LABELS) +
  scale_fill_manual(
    values = METHOD_GROUP_FILLS,
    breaks = c("Reference", "Proposed"),
    labels = c("Reference methods", "Proposed methods")
  ) +
  scale_colour_manual(
    values = METHOD_GROUP_POINTS,
    breaks = c("Reference", "Proposed"),
    labels = c("Reference methods", "Proposed methods")
  ) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  coord_cartesian(ylim = c(-0.05, 1.02)) +
  labs(x = NULL, y = "Subspace similarity (PC10)", fill = NULL) +
  guides(colour = "none") +
  theme_journal(6.2) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.text = element_text(size = 6.2),
    legend.key.width = unit(5.0, "mm"),
    legend.spacing.x = unit(1.2, "mm")
  )
save_pub_r(fig07, fig_path("Real_data", "Figure07_realdata_subspace_pc10"), 183, 78)

if (figure_only == "Figure06_07") {
  cat("Regenerated Figures 6 and 7 from Source Data.\n")
  quit(save = "no", status = 0L)
}

sensitivity <- prepare_methods(read_source("FigureS07_realdata_pc10_vs_pc20_gmm_reproducibility_source_data.csv"))
fig_s7 <- ggplot(sensitivity, aes(y = Method)) +
  geom_segment(aes(x = PC10, xend = PC20, yend = Method),
               colour = "#A7A7A7", linewidth = 0.45) +
  geom_point(aes(x = PC10, shape = "PC10"), colour = "#0072B2", size = 1.8) +
  geom_point(aes(x = PC20, shape = "PC20"), colour = "#D55E00", size = 1.8) +
  facet_wrap(~DatasetName, nrow = 1) +
  scale_y_discrete(limits = rev(METHOD_ORDER), labels = METHOD_LABELS) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  scale_shape_manual(values = c(PC10 = 16, PC20 = 17), name = NULL) +
  labs(x = "GMM FULL-vs-subset ARI", y = NULL) +
  theme_journal(6.2) +
  theme(legend.position = "top")
save_pub_r(
  fig_s7,
  fig_path("Appendix", "FigureS07_realdata_pc10_vs_pc20_gmm_reproducibility"),
  183, 78
)

cat("Regenerated real-data figures from curated Source Data.\n")
