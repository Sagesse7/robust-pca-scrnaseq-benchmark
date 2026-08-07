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

source_dir <- file.path(release_root, "Source_Data")
figure_dir <- file.path(release_root, "figures", "Appendix")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

replicate_file <- file.path(
  source_dir, "FigureS04_simulation2_pcwise_pc10_replicate_source_data.csv"
)
rank_file <- file.path(
  source_dir, "FigureS05_simulation2_pcwise_rank_source_data.csv"
)
if (!file.exists(replicate_file) || !file.exists(rank_file)) {
  stop("Simulation 2 PC-wise Source Data are missing.")
}

replicate_dt <- fread(replicate_file)
rank_dt <- fread(rank_file)
if ("PerturbationCondition" %in% names(rank_dt)) {
  setnames(rank_dt, "PerturbationCondition", "Noise")
}
rank_dt[Noise == "Affected-cell fraction", Noise := "Shifted-cell fraction"]
rank_dt[Noise == "Shift factor", Noise := "Shift magnitude"]

sweeps <- list(
  dropout = list(
    param_ids = 2:6,
    level_column = "dropout_mid",
    levels = as.character(c(-1, 0, 1, 2, 3)),
    labels = c("-1", "0", "1", "2", "3"),
    source_name = "Dropout midpoint",
    title = "Dropout",
    xlab = "Dropout midpoint"
  ),
  doublet = list(
    param_ids = 7:11,
    level_column = "doublet_frac",
    levels = as.character(c(0.02, 0.05, 0.08, 0.10, 0.15)),
    labels = c("2%", "5%", "8%", "10%", "15%"),
    source_name = "Doublet fraction",
    title = "Synthetic doublets",
    xlab = "Doublet fraction"
  ),
  shift_fraction = list(
    param_ids = 12:16,
    level_column = "cell_frac",
    levels = as.character(c(0.02, 0.04, 0.06, 0.08, 0.10)),
    labels = c("2%", "4%", "6%", "8%", "10%"),
    source_name = "Shifted-cell fraction",
    title = "Gene-subset mean shift",
    xlab = "Shifted-cell fraction"
  ),
  shift_magnitude = list(
    param_ids = 17:21,
    level_column = "meanlog",
    levels = as.character(c(1.5, 2, 3, 4, 5)),
    labels = c("1.5", "2", "3", "4", "5"),
    source_name = "Shift magnitude",
    title = "Gene-subset mean shift",
    xlab = "Shift magnitude"
  )
)

distribution_panel <- function(spec, key) {
  z <- copy(replicate_dt[ParamID %in% spec$param_ids & is.finite(PCwisePC10)])
  z[, Method := factor(clean_method(Method), levels = METHOD_ORDER)]
  z[, Level := factor(
    as.character(get(spec$level_column)),
    levels = spec$levels,
    labels = spec$labels
  )]

  reference <- z[0]
  if (key == "dropout") {
    reference <- z[as.character(Level) == "0", .(
      PCwisePC10 = mean(PCwisePC10)
    ), by = Method]
    reference[, Level := factor("0", levels = spec$labels)]
    z <- z[as.character(Level) != "0"]
  }

  axis_labels <- spec$labels
  axis_title <- spec$xlab
  if (key %in% c("doublet", "shift_fraction")) {
    axis_labels <- sub("%$", "", axis_labels)
    axis_title <- paste0(axis_title, " (%)")
  }

  p <- ggplot(z, aes(Level, PCwisePC10)) +
    geom_boxplot(
      width = 0.62, outlier.shape = NA, fill = "#D7E0EA",
      colour = "#333333", linewidth = 0.28
    ) +
    geom_point(
      position = position_jitter(
        width = 0.10, height = 0,
        seed = c(dropout = 451L, doublet = 452L,
                 shift_fraction = 453L, shift_magnitude = 454L)[[key]]
      ),
      size = 0.28, alpha = 0.22, colour = "#222222"
    ) +
    facet_wrap(~Method, ncol = 5, labeller = as_labeller(METHOD_LABELS)) +
    coord_cartesian(ylim = c(0, 1.02)) +
    scale_x_discrete(limits = spec$labels, labels = axis_labels, drop = FALSE) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    labs(x = axis_title, y = "PC-wise similarity (PC1-PC10)", title = spec$title) +
    theme_journal(5.7) +
    theme(
      legend.position = "none",
      strip.text = element_text(size = 5.5),
      axis.text.x = element_text(size = 4.8),
      axis.text.y = element_text(size = 5.0),
      axis.title = element_text(size = 5.4),
      plot.title = element_text(size = 6.2, face = "bold"),
      panel.spacing = grid::unit(1.2, "mm")
    )

  if (nrow(reference)) {
    reference_position <- match("0", spec$labels)
    p <- p +
      annotate(
        "rect", xmin = reference_position - 0.5,
        xmax = reference_position + 0.5,
        ymin = -Inf, ymax = Inf, fill = "#E2E2E2", alpha = 0.55
      ) +
      geom_point(
        data = reference, aes(Level, PCwisePC10), inherit.aes = FALSE,
        shape = 21, fill = "white", colour = "#555555",
        size = 1.05, stroke = 0.3
      )
  }
  p
}

distribution_panels <- Map(distribution_panel, sweeps, names(sweeps))
figure_s4 <- wrap_plots(distribution_panels, ncol = 2) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))
save_pub_r(
  figure_s4,
  file.path(figure_dir, "FigureS04_simulation2_pcwise_pc10"),
  183, 156
)

rank_panel <- function(spec) {
  z <- copy(rank_dt[Noise == spec$source_name])
  z[, `:=`(
    Method = factor(clean_method(Method), levels = METHOD_ORDER),
    Level = factor(as.character(Level), levels = spec$labels),
    Label = ifelse(is.na(Rank), "NR", sprintf("%.0f", Rank)),
    TextColour = ifelse(!is.na(Rank) & Rank >= 6, "white", "#151515")
  )]

  ggplot(z, aes(Level, Method, fill = Rank)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(aes(label = Label, colour = TextColour), size = 1.8) +
    scale_colour_identity() +
    scale_y_discrete(limits = rev(METHOD_ORDER), labels = METHOD_LABELS) +
    scale_fill_viridis_c(
      option = "E", direction = -1, limits = c(1, 10),
      breaks = c(1, 5, 10), na.value = "#E2E2E2", name = "Rank"
    ) +
    labs(title = spec$title, x = spec$xlab, y = NULL) +
    theme_heatmap(6.2) +
    theme(
      axis.text.x = element_text(size = 5.3),
      axis.text.y = element_text(size = 5.3),
      plot.title = element_text(size = 6.5, face = "bold")
    )
}

rank_panels <- lapply(sweeps, rank_panel)
figure_s5 <- wrap_plots(rank_panels, ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "right", plot.tag = element_text(size = 8, face = "bold"))
save_pub_r(
  figure_s5,
  file.path(figure_dir, "FigureS05_simulation2_pcwise_rank"),
  183, 112
)

cat("Regenerated Simulation 2 PC-wise figures from curated Source Data.\n")
