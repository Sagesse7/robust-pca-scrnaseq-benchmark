#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(svglite)
})

script_file <- sub(
  "^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
)
script_dir <- dirname(normalizePath(script_file))
release_root <- dirname(dirname(script_dir))
source(file.path(release_root, "scripts", "utils", "figure_style.R"))

parse_cli <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) next
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    out[[parts[1]]] <- paste(parts[-1], collapse = "=")
  }
  out
}

cli <- parse_cli(commandArgs(trailingOnly = TRUE))
source_file <- if (is.null(cli[["source-file"]])) {
  file.path(release_root, "Source_Data", "Figure08_runtime_scaling_source_data.csv")
} else cli[["source-file"]]
output_dir <- if (is.null(cli[["output-dir"]])) {
  file.path(release_root, "figures", "Runtime")
} else cli[["output-dir"]]
expected_repeats <- if (is.null(cli[["expected-repeats"]])) {
  5L
} else {
  as.integer(cli[["expected-repeats"]])
}

if (!file.exists(source_file)) stop("Runtime Source Data not found: ", source_file)
runtime <- fread(source_file)
required <- c(
  "MethodLabel", "Cell", "Gene", "CellType", "N", "MedianRuntimeSec",
  "Q25RuntimeSec", "Q75RuntimeSec", "Calibration", "CodeVersion",
  "FigurePanel"
)
missing_columns <- setdiff(required, names(runtime))
if (length(missing_columns)) {
  stop("Runtime summary is missing columns: ", paste(missing_columns, collapse = ", "))
}
if (expected_repeats < 1L || any(runtime$N != expected_repeats)) {
  stop("Every plotted setting/method must contain the requested number of repeats.")
}
if (uniqueN(runtime$CodeVersion) != 1L) stop("Mixed runtime code versions detected.")

runtime[, Method := clean_method(MethodLabel)]
runtime[, Method := factor(Method, levels = METHOD_ORDER)]

runtime_panel <- function(z, xvar, breaks, xlab, panel_title) {
  ggplot(
    z,
    aes(
      x = get(xvar), y = MedianRuntimeSec, colour = Method,
      shape = Method, linetype = Method, group = Method
    )
  ) +
    geom_errorbar(
      aes(ymin = Q25RuntimeSec, ymax = Q75RuntimeSec),
      width = 35, linewidth = 0.30, alpha = 0.55
    ) +
    geom_line(linewidth = 0.62) +
    geom_point(size = 1.35, stroke = 0.32) +
    scale_method_colour() +
    scale_method_shape() +
    scale_method_linetype() +
    scale_x_continuous(breaks = breaks) +
    scale_y_log10(
      breaks = c(0.1, 1, 10, 100, 1000, 10000),
      labels = label_log()
    ) +
    labs(
      x = xlab, y = "Runtime (s)", title = panel_title,
      colour = NULL, shape = NULL, linetype = NULL
    ) +
    theme_journal(6.2) +
    theme(
      legend.position = "bottom",
      legend.key.width = unit(5.2, "mm"),
      plot.title = element_text(size = 7, face = "bold", hjust = 0)
    ) +
    guides(
      colour = guide_legend(nrow = 2, byrow = TRUE),
      shape = guide_legend(nrow = 2, byrow = TRUE),
      linetype = guide_legend(nrow = 2, byrow = TRUE)
    )
}

gene_data <- runtime[FigurePanel == "A_gene_scaling"]
cell_data <- runtime[FigurePanel == "B_cell_scaling"]
if (nrow(gene_data) != 40L || nrow(cell_data) != 40L) {
  stop("Expected four scaling settings and ten methods in each runtime panel.")
}

gene_panel <- runtime_panel(
  gene_data, "Gene", c(500, 1000, 2000, 3000), "Genes, d",
  "Runtime across gene dimensions (n = 2,000)"
)
cell_panel <- runtime_panel(
  cell_data, "Cell", c(500, 1000, 2000, 3000), "Cells, n",
  "Runtime across cell numbers (d = 2,000)"
) + labs(y = NULL)

combined <- (gene_panel | cell_panel) +
  plot_layout(guides = "collect", widths = c(1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(
    legend.position = "bottom",
    legend.key.width = unit(5.4, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    plot.tag = element_text(size = 8.5, face = "bold")
  )

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
stem <- file.path(output_dir, "Figure08_runtime_scaling")
save_pub_r(combined, stem, 183, 86)
cat("Regenerated matched-design runtime figure from curated Source Data at:",
    normalizePath(output_dir), "\n")
