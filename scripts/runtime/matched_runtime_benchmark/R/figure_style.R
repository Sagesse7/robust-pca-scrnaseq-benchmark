METHOD_ORDER <- c(
  "PCA", "Grid", "Hubert", "PCP", "K's tau",
  "Winsor", "Quad", "Ball", "Shell", "LR"
)

# Internal computation IDs remain Grid/Hubert/Tau-compatible.  All manuscript
# figures and publication Source Data use the labels below, in this order.
PUBLICATION_METHOD_ORDER <- c(
  "PCA", "PcaGrid", "PcaHubert", "PCP", "K's tau",
  "Winsor", "Quad", "Ball", "Shell", "LR"
)

METHOD_LABELS <- c(
  "PCA" = "PCA",
  "Grid" = "PcaGrid",
  "Hubert" = "PcaHubert",
  "PCP" = "PCP",
  "K's tau" = "K's tau",
  "Winsor" = "Winsor",
  "Quad" = "Quad",
  "Ball" = "Ball",
  "Shell" = "Shell",
  "LR" = "LR"
)

METHOD_COLORS <- c(
  "PCA" = "#222222",
  "Grid" = "#777777",
  "Hubert" = "#3B6FB6",
  "PCP" = "#6F90C8",
  "K's tau" = "#2A9D8F",
  "Winsor" = "#E69F00",
  "Quad" = "#D55E00",
  "Ball" = "#CC79A7",
  "Shell" = "#A75D91",
  "LR" = "#7A3E65"
)

METHOD_SHAPES <- c(
  "PCA" = 16, "Grid" = 4, "Hubert" = 17, "PCP" = 15,
  "K's tau" = 18, "Winsor" = 1, "Quad" = 2, "Ball" = 0,
  "Shell" = 5, "LR" = 8
)

METHOD_LINETYPES <- c(
  "PCA" = "solid", "Grid" = "22", "Hubert" = "solid",
  "PCP" = "longdash", "K's tau" = "dotdash", "Winsor" = "solid",
  "Quad" = "longdash", "Ball" = "solid", "Shell" = "dotdash",
  "LR" = "dashed"
)

PROPOSED_METHODS <- c("Winsor", "Quad", "Ball", "Shell", "LR")

clean_method <- function(x) {
  x <- as.character(x)
  x[x == "PcaGrid"] <- "Grid"
  x[x == "PcaHubert"] <- "Hubert"
  x[x == "Tau"] <- "K's tau"
  factor(x, levels = METHOD_ORDER)
}

display_method <- function(x) {
  unname(METHOD_LABELS[as.character(x)])
}

publication_method <- function(x) {
  internal <- clean_method(x)
  unname(METHOD_LABELS[as.character(internal)])
}

publication_source_data <- function(x, group_cols = character()) {
  z <- data.table::as.data.table(data.table::copy(x))
  if (!"Method" %in% names(z)) return(z)
  z[, Method := publication_method(Method)]
  z[, publication_method_order__ := match(Method, PUBLICATION_METHOD_ORDER)]
  order_cols <- c(intersect(group_cols, names(z)), "publication_method_order__")
  data.table::setorderv(z, order_cols, na.last = TRUE)
  z[, publication_method_order__ := NULL]
  z[]
}

method_group <- function(x) {
  ifelse(as.character(x) %in% PROPOSED_METHODS, "Proposed", "Reference")
}

theme_journal <- function(base_size = 6.5, base_family = "Arial") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = "#252525"),
      axis.ticks = ggplot2::element_line(linewidth = 0.35, colour = "#252525"),
      axis.ticks.length = grid::unit(1.4, "mm"),
      axis.title = ggplot2::element_text(size = base_size, colour = "#202020"),
      axis.text = ggplot2::element_text(size = base_size - 0.5, colour = "#202020"),
      legend.title = ggplot2::element_text(size = base_size - 0.2),
      legend.text = ggplot2::element_text(size = base_size - 0.5),
      legend.key.height = grid::unit(3.2, "mm"),
      legend.key.width = grid::unit(4.2, "mm"),
      legend.background = ggplot2::element_blank(),
      legend.box.background = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "#F1F1F1", colour = NA),
      strip.text = ggplot2::element_text(size = base_size, face = "bold", colour = "#202020"),
      plot.title = ggplot2::element_text(size = base_size + 0.5, face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(size = base_size, hjust = 0),
      plot.tag = ggplot2::element_text(size = 8, face = "bold"),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(4, 5, 4, 5)
    )
}

theme_heatmap <- function(base_size = 6.5, base_family = "Arial") {
  theme_journal(base_size, base_family) +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, vjust = 1),
      plot.margin = ggplot2::margin(3, 4, 3, 4)
    )
}

save_pub_r <- function(plot, stem, width_mm, height_mm, bg = "white") {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4

  svg_file <- paste0(stem, ".svg")
  pdf_file <- paste0(stem, ".pdf")

  svglite::svglite(
    svg_file, width = width_in, height = height_in,
    bg = bg, system_fonts = list(sans = "Arial")
  )
  print(plot)
  grDevices::dev.off()

  # This R installation lacks the XQuartz libraries required by cairo_pdf().
  # Convert the canonical editable SVG with librsvg instead; artwork stays vector.
  rsvg::rsvg_pdf(svg_file, pdf_file)

  invisible(c(pdf = pdf_file, svg = svg_file))
}

scale_method_colour <- function(...) {
  ggplot2::scale_colour_manual(
    values = METHOD_COLORS,
    breaks = METHOD_ORDER,
    labels = METHOD_LABELS[METHOD_ORDER],
    drop = FALSE,
    ...
  )
}

scale_method_fill <- function(...) {
  ggplot2::scale_fill_manual(
    values = METHOD_COLORS,
    breaks = METHOD_ORDER,
    labels = METHOD_LABELS[METHOD_ORDER],
    drop = FALSE,
    ...
  )
}

scale_method_shape <- function(...) {
  ggplot2::scale_shape_manual(
    values = METHOD_SHAPES,
    breaks = METHOD_ORDER,
    labels = METHOD_LABELS[METHOD_ORDER],
    drop = FALSE,
    ...
  )
}

scale_method_linetype <- function(...) {
  ggplot2::scale_linetype_manual(
    values = METHOD_LINETYPES,
    breaks = METHOD_ORDER,
    labels = METHOD_LABELS[METHOD_ORDER],
    drop = FALSE,
    ...
  )
}
