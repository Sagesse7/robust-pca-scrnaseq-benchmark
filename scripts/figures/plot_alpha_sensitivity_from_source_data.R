suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
  library(rsvg)
})

script_file <- sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)[1])
script_dir <- dirname(normalizePath(script_file))
release_root <- dirname(dirname(script_dir))
source(file.path(release_root, "scripts", "utils", "figure_style.R"))

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) args[1] else file.path(release_root, "figures", "Appendix")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

z <- fread(file.path(release_root, "Source_Data", "FigureS08_alpha_sensitivity_source_data.csv"))
z[, Method := factor(Method, levels = PROPOSED_METHODS)]
z[, AlphaLabel := factor(
  sprintf("%.2f", Alpha), levels = c("0.10", "0.20"),
  labels = c(expression(alpha == 0.10), expression(alpha == 0.20))
)]

alpha_colours <- c("0.10" = "#4C78A8", "0.20" = "#E28E2C")

delta_panel <- function(data, outcomes, y_limits, y_breaks, title) {
  pdat <- copy(data[Outcome %chin% outcomes])
  pdat[, Outcome := factor(Outcome, levels = outcomes)]
  ggplot(pdat, aes(Method, Delta, fill = sprintf("%.2f", Alpha))) +
    geom_hline(yintercept = 0, linewidth = 0.32, colour = "#4A4A4A") +
    geom_boxplot(
      position = position_dodge2(width = 0.70, preserve = "single"),
      width = 0.58, outlier.shape = NA, linewidth = 0.28, alpha = 0.82
    ) +
    geom_point(
      aes(colour = sprintf("%.2f", Alpha)),
      position = position_jitterdodge(
        jitter.width = 0.10, jitter.height = 0, dodge.width = 0.70, seed = 711L
      ),
      size = 0.42, alpha = 0.38, stroke = 0
    ) +
    facet_wrap(~Outcome, nrow = 1) +
    scale_fill_manual(values = alpha_colours, name = expression(alpha)) +
    scale_colour_manual(values = alpha_colours, guide = "none") +
    scale_y_continuous(limits = y_limits, breaks = y_breaks) +
    labs(
      title = title, x = NULL,
      y = expression(Delta~"versus"~alpha == 0.05)
    ) +
    theme_journal(6.3) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "top",
      legend.justification = "left",
      strip.text = element_text(size = 6.1, face = "bold"),
      plot.title = element_text(size = 6.8, face = "bold")
    )
}

p_a <- delta_panel(z, "GMM ARI", c(-0.075, 0.34), c(-0.05, 0, 0.1, 0.2, 0.3),
                   "GMM clustering sensitivity")
p_b <- delta_panel(z, "GMM NMI", c(-0.06, 0.27), c(-0.05, 0, 0.1, 0.2),
                   "GMM clustering sensitivity")
p_c <- delta_panel(z, c("k-means ARI", "Louvain ARI"), c(-0.055, 0.055),
                   c(-0.05, 0, 0.05), "Alternative clustering algorithms")
p_d <- delta_panel(z, c("Subspace PC10", "Subspace PC20"), c(-0.055, 0.055),
                   c(-0.05, 0, 0.05), "Subspace sensitivity")

fig <- (p_a | p_b) / (p_c | p_d) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(
    legend.position = "top",
    plot.tag = element_text(size = 8, face = "bold")
  )

stem <- file.path(output_dir, "FigureS08_alpha_sensitivity")
save_pub_r(fig, stem, 183, 132)
ragg::agg_png(
  paste0(stem, ".png"), width = 183, height = 132, units = "mm", res = 300,
  background = "white"
)
print(fig)
dev.off()

cat("Wrote Figure S8 to", output_dir, "\n")
