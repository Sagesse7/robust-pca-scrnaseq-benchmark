suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(svglite)
})

script_file <- sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)[1])
script_dir <- dirname(normalizePath(script_file))
release_root <- dirname(dirname(script_dir))
source(file.path(release_root, "scripts", "utils", "figure_style.R"))

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1]]) else default
}

source_file <- arg_value(
  "source-file",
  file.path(release_root, "Source_Data", "simulation_figures_input_alpha005.csv")
)
out_stem <- arg_value(
  "out-stem",
  file.path(
    release_root, "figures", "Simulation_1",
    "Figure03_simulation1_gmm_pc10_clustering_ari_nmi"
  )
)
source_out <- arg_value(
  "source-out",
  file.path(
    release_root, "Source_Data",
    "Figure03_simulation1_gmm_pc10_clustering_ari_nmi_source_data.csv"
  )
)

if (!file.exists(source_file)) stop("Source data not found: ", source_file)

sim1 <- fread(source_file)
sim1 <- sim1[is.na(Alpha) | abs(Alpha - 0.05) < 1e-12]
sim1[, Method := clean_method(Method)]

noise_specs <- list(
  dropout = list(
    param_ids = 2:6,
    levels = as.character(c(-1, 0, 1, 2, 3)),
    labels = c("-1", "0", "1", "2", "3"),
    title = "Dropout midpoint",
    column = "dropout_mid"
  ),
  doublet = list(
    param_ids = 7:11,
    levels = as.character(c(0.02, 0.05, 0.08, 0.10, 0.15)),
    labels = c("2%", "5%", "8%", "10%", "15%"),
    title = "Doublet fraction",
    column = "doublet_frac"
  ),
  shift_fraction = list(
    param_ids = 12:16,
    levels = as.character(c(0.02, 0.04, 0.06, 0.08, 0.10)),
    labels = c("2%", "4%", "6%", "8%", "10%"),
    title = "Affected-cell fraction",
    column = "cell_frac"
  ),
  shift_factor = list(
    param_ids = 17:21,
    levels = as.character(c(1.5, 2, 3, 4, 5)),
    labels = c("1.5", "2", "3", "4", "5"),
    title = "Mean-log shift",
    column = "meanlog"
  )
)

metric_specs <- list(
  ARI = "ARI_gmm_pc10",
  NMI = "NMI_gmm_pc10"
)

plot_dt <- rbindlist(lapply(names(noise_specs), function(noise_key) {
  spec <- noise_specs[[noise_key]]
  z <- copy(sim1[ParamID %in% spec$param_ids])
  z[, LevelRaw := as.character(get(spec$column))]
  level_map <- data.table(
    LevelRaw = spec$levels,
    LevelKey = paste(noise_key, seq_along(spec$levels), sep = "__"),
    LevelLabel = spec$labels
  )
  z <- merge(z, level_map, by = "LevelRaw", all.x = TRUE, sort = FALSE)
  z[, Noise := spec$title]
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

noise_order <- vapply(noise_specs, `[[`, character(1), "title")
level_order <- unlist(lapply(names(noise_specs), function(noise_key) {
  paste(noise_key, seq_along(noise_specs[[noise_key]]$levels), sep = "__")
}), use.names = FALSE)
level_labels <- unlist(lapply(names(noise_specs), function(noise_key) {
  setNames(noise_specs[[noise_key]]$labels, paste(noise_key, seq_along(noise_specs[[noise_key]]$levels), sep = "__"))
}))
plot_dt[, Noise := factor(Noise, levels = noise_order)]
plot_dt[, Metric := factor(Metric, levels = c("ARI", "NMI"))]
plot_dt[, LevelKey := factor(LevelKey, levels = level_order)]
plot_dt[, Method := factor(as.character(Method), levels = METHOD_ORDER)]
plot_dt[, Label := sprintf("%.2f", Mean)]
plot_dt[, TextColour := ifelse(Mean < 0.43, "white", "#151515")]

source_dt <- copy(plot_dt)
source_dt[, Method := publication_method(Method)]
setnames(source_dt, "LevelLabel", "Level")
source_dt <- source_dt[, .(Metric, Noise, Level, Method, Mean, SD, N)]
fwrite(source_dt, source_out)

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
    option = "E", limits = c(0, 1), oob = squish,
    breaks = c(0, 0.5, 1), name = "Score"
  ) +
  facet_grid(
    Metric ~ Noise, scales = "free_x", space = "free_x",
    switch = "y"
  ) +
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

save_pub_r(fig, out_stem, width_mm = 183, height_mm = 128)

png_file <- paste0(out_stem, ".png")
# Render from the canonical SVG. This avoids Poppler filename-format warnings
# and does not introduce an additional PDF-rendering dependency.
rsvg::rsvg_png(paste0(out_stem, ".svg"), png_file, width = 1585)

message("Wrote: ", paste0(out_stem, ".pdf"))
message("Wrote: ", paste0(out_stem, ".svg"))
message("Wrote: ", png_file)
message("Wrote: ", source_out)
