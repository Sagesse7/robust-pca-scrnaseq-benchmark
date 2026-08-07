suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(rsvg)
})

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

sim1 <- read_source("simulation_figures_input_alpha005.csv")
sim2 <- read_source("subspace_figures_input_alpha005.csv")
sim1 <- sim1[is.na(Alpha) | abs(Alpha - 0.05) < 1e-12]
sim2 <- sim2[is.na(Alpha) | abs(Alpha - 0.05) < 1e-12]
sim1[, Method := clean_method(Method)]
sim2[, Method := clean_method(Method)]

make_sweeps <- function(dt) {
  list(
    dropout = list(
      data = copy(dt[ParamID %in% 2:6])[, level := as.character(dropout_mid)],
      levels = as.character(c(-1, 0, 1, 2, 3)),
      labels = c("-1", "0", "1", "2", "3"),
      title = "Dropout", xlab = "Dropout midpoint"
    ),
    doublet = list(
      data = copy(dt[ParamID %in% 7:11])[, level := as.character(doublet_frac)],
      levels = as.character(c(0.02, 0.05, 0.08, 0.10, 0.15)),
      labels = c("2%", "5%", "8%", "10%", "15%"),
      title = "Synthetic doublets", xlab = "Doublet fraction"
    ),
    shift_fraction = list(
      data = copy(dt[ParamID %in% 12:16])[, level := as.character(cell_frac)],
      levels = as.character(c(0.02, 0.04, 0.06, 0.08, 0.10)),
      labels = c("2%", "4%", "6%", "8%", "10%"),
      title = "Gene-subset mean shift", xlab = "Shifted-cell fraction"
    ),
    shift_factor = list(
      data = copy(dt[ParamID %in% 17:21])[, level := as.character(meanlog)],
      levels = as.character(c(1.5, 2, 3, 4, 5)),
      labels = c("1.5", "2", "3", "4", "5"),
      title = "Gene-subset mean shift", xlab = "Shift magnitude"
    )
  )
}

sim1_sweeps <- make_sweeps(sim1)
sim2_sweeps <- make_sweeps(sim2)

build_noise_summary <- function(metric_specs) {
  plot_dt <- rbindlist(lapply(names(sim1_sweeps), function(noise_key) {
    spec <- sim1_sweeps[[noise_key]]
    z <- copy(spec$data)
    z[, LevelRaw := level]
    level_map <- data.table(
      LevelRaw = spec$levels,
      LevelKey = paste(noise_key, seq_along(spec$levels), sep = "__"),
      LevelLabel = spec$labels
    )
    z <- merge(z, level_map, by = "LevelRaw", all.x = TRUE, sort = FALSE)
    z[, Noise := spec$xlab]
    rbindlist(lapply(names(metric_specs), function(metric_name) {
      metric_col <- metric_specs[[metric_name]]
      zz <- z[, .(
        Mean = mean(get(metric_col), na.rm = TRUE),
        SD = sd(get(metric_col), na.rm = TRUE),
        N = sum(is.finite(get(metric_col)))
      ), by = .(Noise, LevelKey, LevelLabel, Method)]
      zz[, Metric := metric_name]
      zz
    }))
  }), use.names = TRUE)

  noise_order <- vapply(sim1_sweeps, `[[`, character(1), "xlab")
  level_order <- unlist(lapply(names(sim1_sweeps), function(noise_key) {
    paste(noise_key, seq_along(sim1_sweeps[[noise_key]]$levels), sep = "__")
  }), use.names = FALSE)
  level_labels <- unlist(lapply(names(sim1_sweeps), function(noise_key) {
    setNames(
      sim1_sweeps[[noise_key]]$labels,
      paste(noise_key, seq_along(sim1_sweeps[[noise_key]]$levels), sep = "__")
    )
  }))

  plot_dt[, `:=`(
    Noise = factor(Noise, levels = noise_order),
    Metric = factor(Metric, levels = names(metric_specs)),
    LevelKey = factor(LevelKey, levels = level_order),
    Method = factor(as.character(Method), levels = METHOD_ORDER)
  )]
  list(data = plot_dt, level_labels = level_labels)
}

save_combined_clustering <- function(metric_specs, figure_stem, source_filename) {
  built <- build_noise_summary(metric_specs)
  plot_dt <- built$data
  level_labels <- built$level_labels
  plot_dt[, `:=`(
    Label = sprintf("%.2f", fifelse(abs(Mean) < 0.005, 0, Mean)),
    TextColour = ifelse(Mean < 0.43, "white", "#151515")
  )]

  source_dt <- copy(plot_dt)
  source_dt[, `:=`(
    MetricOrder = as.integer(Metric),
    NoiseOrder = as.integer(Noise),
    LevelOrder = as.integer(LevelKey),
    MethodOrder = as.integer(Method)
  )]
  source_dt[, Method := publication_method(Method)]
  setorder(source_dt, MetricOrder, NoiseOrder, LevelOrder, MethodOrder)
  setnames(source_dt, "LevelLabel", "Level")
  source_dt <- source_dt[, .(Metric, Noise, Level, Method, Mean, SD, N)]
  setnames(source_dt, "Noise", "PerturbationCondition")
  fwrite(source_dt, file.path(source_data_dir, source_filename))

  fig <- ggplot(plot_dt, aes(LevelKey, Method, fill = Mean)) +
    geom_tile(colour = "white", linewidth = 0.22) +
    geom_text(
      aes(label = Label, colour = TextColour),
      size = 1.55, lineheight = 0.85, show.legend = FALSE
    ) +
    scale_colour_identity() +
    scale_y_discrete(limits = rev(METHOD_ORDER), labels = METHOD_LABELS) +
    scale_x_discrete(labels = level_labels, drop = TRUE) +
    scale_fill_viridis_c(
      option = "E", limits = c(0, 1), oob = scales::squish,
      breaks = c(0, 0.5, 1), name = "Score"
    ) +
    facet_grid(Metric ~ Noise, scales = "free_x", space = "free_x", switch = "y") +
    labs(x = NULL, y = NULL) +
    theme_heatmap(6.2) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 6.2),
      legend.text = element_text(size = 5.8),
      strip.placement = "outside",
      strip.background = element_rect(fill = "#F1F1F1", colour = NA),
      strip.text.x = element_text(size = 6.4, face = "bold"),
      strip.text.y.left = element_text(size = 6.4, face = "bold", angle = 0),
      axis.text.x = element_text(size = 5.1, angle = 35, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 5.6),
      panel.spacing.x = grid::unit(1.2, "mm"),
      panel.spacing.y = grid::unit(1.3, "mm"),
      plot.margin = margin(3, 4, 3, 3)
    )

  save_pub_r(fig, fig_path("Appendix", figure_stem), 183, 128)
}

ari_built <- build_noise_summary(c(ARI = "ARI_gmm_pc10"))
g_built <- build_noise_summary(c(`Selected G` = "k_gmm_pc10"))
level_labels <- ari_built$level_labels
ari_dt <- copy(ari_built$data)
g_dt <- copy(g_built$data)
setnames(ari_dt, c("Mean", "SD", "N"), c("MeanARI", "SDARI", "NARI"))
setnames(g_dt, c("Mean", "SD", "N"), c("MeanG", "SDG", "NG"))
s1_dt <- merge(
  ari_dt[, .(Noise, LevelKey, LevelLabel, Method, MeanARI, SDARI, NARI)],
  g_dt[, .(Noise, LevelKey, Method, MeanG, SDG, NG)],
  by = c("Noise", "LevelKey", "Method"), all = TRUE, sort = FALSE
)
s1_dt[, `:=`(
  Label = sprintf("%.2f\nG=%.1f", fifelse(abs(MeanARI) < 0.005, 0, MeanARI), MeanG),
  TextColour = ifelse(MeanARI < 0.43, "white", "#151515")
)]

s1_source <- copy(s1_dt)
s1_source[, `:=`(
  NoiseOrder = as.integer(Noise),
  LevelOrder = as.integer(LevelKey),
  MethodOrder = as.integer(Method)
)]
s1_source[, Method := publication_method(Method)]
setorder(s1_source, NoiseOrder, LevelOrder, MethodOrder)
setnames(s1_source, "LevelLabel", "Level")
s1_source <- s1_source[, .(
  Noise, Level, Method, MeanARI, SDARI, MeanG, SDG,
  N = pmin(NARI, NG)
)]
setnames(s1_source, "Noise", "PerturbationCondition")
fwrite(
  s1_source,
  file.path(source_data_dir, "FigureS01_simulation1_gmm_pc10_ari_with_k_source_data.csv")
)

fig_s1 <- ggplot(s1_dt, aes(LevelKey, Method, fill = MeanARI)) +
  geom_tile(colour = "white", linewidth = 0.22) +
  geom_text(
    aes(label = Label, colour = TextColour),
    size = 1.42, lineheight = 0.86, show.legend = FALSE
  ) +
  scale_colour_identity() +
  scale_y_discrete(limits = rev(METHOD_ORDER), labels = METHOD_LABELS) +
  scale_x_discrete(labels = level_labels, drop = TRUE) +
  scale_fill_viridis_c(
    option = "E", limits = c(0, 1), oob = scales::squish,
    breaks = c(0, 0.5, 1), name = "ARI"
  ) +
  facet_grid(. ~ Noise, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = NULL) +
  theme_heatmap(6.2) +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "#F1F1F1", colour = NA),
    strip.text.x = element_text(size = 6.4, face = "bold"),
    axis.text.x = element_text(size = 5.1, angle = 35, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 5.6),
    panel.spacing.x = grid::unit(1.2, "mm"),
    plot.margin = margin(3, 4, 3, 3)
  )
save_pub_r(
  fig_s1,
  fig_path("Appendix", "FigureS01_simulation1_gmm_pc10_ari_with_k"),
  183, 82
)

if ("--figure-s1-only" %in% commandArgs(trailingOnly = TRUE)) {
  cat("Regenerated Figure S1 from curated Source Data.\n")
  quit(save = "no", status = 0)
}

save_combined_clustering(
  c(ARI = "ARI_kmeans_pc10", NMI = "NMI_kmeans_pc10"),
  "FigureS02_simulation1_kmeans_pc10_ari_nmi",
  "FigureS02_simulation1_kmeans_pc10_ari_nmi_source_data.csv"
)
save_combined_clustering(
  c(ARI = "ARI_louvain_pc10", NMI = "NMI_louvain_pc10"),
  "FigureS03_simulation1_louvain_pc10_ari_nmi",
  "FigureS03_simulation1_louvain_pc10_ari_nmi_source_data.csv"
)

distribution_panel <- function(spec, key) {
  z <- copy(spec$data[is.finite(SimilarityPC10)])
  z[, Method := factor(as.character(Method), levels = METHOD_ORDER)]
  z[, Level := factor(level, levels = spec$levels, labels = spec$labels)]
  reference <- z[0]
  if (identical(key, "dropout")) {
    reference <- z[level == "0", .(SimilarityPC10 = mean(SimilarityPC10)), by = Method]
    reference[, Level := factor("0", levels = spec$labels)]
    z <- z[level != "0"]
  }
  axis_labels <- spec$labels
  axis_title <- spec$xlab
  if (key %in% c("doublet", "shift_fraction")) {
    axis_labels <- sub("%$", "", spec$labels)
    axis_title <- paste0(spec$xlab, " (%)")
  }
  jitter_seed <- c(dropout = 401L, doublet = 402L, shift_fraction = 403L, shift_factor = 404L)[[key]]

  p <- ggplot(z, aes(Level, SimilarityPC10)) +
    geom_boxplot(
      width = 0.62, outlier.shape = NA, fill = "#D7E0EA",
      colour = "#333333", linewidth = 0.28
    ) +
    geom_point(
      position = position_jitter(width = 0.10, height = 0, seed = jitter_seed),
      size = 0.28, alpha = 0.22, colour = "#222222"
    ) +
    facet_wrap(~Method, ncol = 5, labeller = as_labeller(METHOD_LABELS)) +
    coord_cartesian(ylim = c(0, 1.02)) +
    scale_x_discrete(limits = spec$labels, labels = axis_labels, drop = FALSE) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    labs(title = spec$title, x = axis_title, y = "Subspace similarity (PC10)") +
    theme_journal(5.2) +
    theme(
      strip.text = element_text(size = 5.2, face = "bold"),
      axis.text = element_text(size = 4.8),
      axis.title = element_text(size = 5.2),
      plot.title = element_text(size = 6.2, face = "bold"),
      panel.spacing = unit(1.2, "mm")
    )

  if (nrow(reference)) {
    reference_position <- match("0", spec$labels)
    p <- p +
      annotate(
        "rect", xmin = reference_position - 0.5, xmax = reference_position + 0.5,
        ymin = -Inf, ymax = Inf, fill = "#E2E2E2", alpha = 0.55
      ) +
      geom_point(
        data = reference, aes(Level, SimilarityPC10), inherit.aes = FALSE,
        shape = 21, fill = "white", colour = "#555555", size = 1.05, stroke = 0.3
      )
  }
  p
}

dist_panels <- Map(distribution_panel, sim2_sweeps, names(sim2_sweeps))
fig04 <- wrap_plots(dist_panels, ncol = 2) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))
save_pub_r(fig04, fig_path("Simulation_2", "Figure04_simulation2_subspace"), 183, 156)

rank_panel <- function(spec, key) {
  z <- copy(spec$data[is.finite(SimilarityPC10)])
  z[, Method := factor(as.character(Method), levels = METHOD_ORDER)]
  z[, Level := factor(level, levels = spec$levels, labels = spec$labels)]
  ranks <- z[, .(Median = median(SimilarityPC10)), by = .(Method, Level)]
  ranks[, Rank := frank(-Median, ties.method = "average"), by = Level]
  if (identical(key, "dropout")) ranks[as.character(Level) == "0", Rank := NA_real_]
  ranks[, Label := ifelse(is.na(Rank), "NR", sprintf("%.0f", Rank))]
  ranks[, TextColour := ifelse(!is.na(Rank) & Rank >= 6, "white", "#151515")]

  ggplot(ranks, aes(Level, Method, fill = Rank)) +
    geom_tile(colour = "white", linewidth = 0.22) +
    geom_text(aes(label = Label, colour = TextColour), size = 1.65, show.legend = FALSE) +
    scale_colour_identity() +
    scale_y_discrete(limits = rev(METHOD_ORDER), labels = METHOD_LABELS) +
    scale_fill_viridis_c(
      option = "E", direction = -1, limits = c(1, 10),
      breaks = c(1, 5, 10), na.value = "#E2E2E2", name = "Rank"
    ) +
    labs(title = spec$title, x = spec$xlab, y = NULL) +
    theme_heatmap(6.2)
}

rank_panels <- Map(rank_panel, sim2_sweeps, names(sim2_sweeps))
fig05 <- wrap_plots(rank_panels, ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "right", plot.tag = element_text(size = 8, face = "bold"))
save_pub_r(fig05, fig_path("Simulation_2", "Figure05_simulation2_rank"), 183, 112)

pc_long <- read_source("FigureS06_simulation1_pc10_vs_pc20_gmm_ari_source_data.csv")
# ParamID 1 and 3 represent the same midpoint-0 background, and ParamID 12
# and 18 represent the shared mean-shift anchor. Retain one representation of
# each scientific condition in the pooled sensitivity summary.
pc_long <- pc_long[ParamID %in% setdiff(2L:21L, 18L)]
stopifnot(
  uniqueN(pc_long$ParamID) == 19L,
  nrow(pc_long) == 3800L
)
pc_long[, Method := clean_method(Method)]
pc_long[, `:=`(
  PC = factor(PC, levels = c("PC10", "PC20")),
  Group = method_group(Method)
)]
fig_s6 <- ggplot(pc_long[is.finite(ARI)], aes(Method, ARI, fill = Group)) +
  geom_boxplot(width = 0.62, outlier.shape = NA, linewidth = 0.32) +
  geom_point(
    position = position_jitter(width = 0.10, height = 0, seed = 408L),
    alpha = 0.17, size = 0.30, colour = "#202020"
  ) +
  facet_wrap(~PC, nrow = 1) +
  scale_x_discrete(labels = METHOD_LABELS) +
  scale_fill_manual(values = c(Reference = "#D4DEE9", Proposed = "#E5C7D5")) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  coord_cartesian(ylim = c(-0.05, 1.02)) +
  labs(x = NULL, y = "Adjusted Rand index", fill = NULL) +
  theme_journal(6.2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
save_pub_r(fig_s6, fig_path("Appendix", "FigureS06_simulation1_pc10_vs_pc20_gmm_ari"), 183, 82)

cat("Regenerated simulation support figures from curated Source Data.\n")
