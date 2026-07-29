suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(ragg)
  library(rsvg)
  library(svglite)
})

script_file <- sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)[1])
script_dir <- dirname(normalizePath(script_file))
release_root <- dirname(dirname(script_dir))
args <- commandArgs(trailingOnly = TRUE)
source_arg <- grep("^--source-data-dir=", args, value = TRUE)
out_arg <- grep("^--out-dir=", args, value = TRUE)
source_data_dir <- if (length(source_arg)) {
  normalizePath(sub("^--source-data-dir=", "", source_arg[[1]]), mustWork = TRUE)
} else {
  file.path(release_root, "Source_Data")
}
out_dir <- if (length(out_arg)) {
  normalizePath(sub("^--out-dir=", "", out_arg[[1]]), mustWork = FALSE)
} else {
  file.path(release_root, "figures", "main_paper")
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_data_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(release_root, "scripts", "utils", "figure_style.R"))

figure2_only <- "--figure2-only" %in% args
figure9_only <- "--figure9-only" %in% args
if (figure2_only && figure9_only) {
  stop("Choose only one of --figure2-only and --figure9-only.")
}

save_drawable <- function(draw, stem, width_mm, height_mm, bg = "white") {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  svg_file <- paste0(stem, ".svg")
  pdf_file <- paste0(stem, ".pdf")
  png_file <- paste0(stem, ".png")

  svglite(svg_file, width = width_in, height = height_in, bg = bg,
          system_fonts = list(sans = "Arial"))
  draw()
  dev.off()
  rsvg_pdf(svg_file, pdf_file)

  agg_png(png_file, width = width_in, height = height_in,
          units = "in", res = 220, background = bg)
  draw()
  dev.off()
  invisible(c(pdf = pdf_file, svg = svg_file, png = png_file))
}

# -----------------------------------------------------------------------------
# Figure 2: weighting functions
# Contract: retain the original Figure 2 scientific message while giving the
# curves more horizontal space and marking the schematic threshold locations.
# -----------------------------------------------------------------------------

q1 <- 0.85
q2 <- 2.15
q3 <- 3.45
q3s <- 4.15
r <- seq(0.01, 5.00, by = 0.001)
pca_reference_weight <- 1 / sqrt(2)

weight_dt <- rbindlist(list(
  data.table(r, Method = "PCA", Weight = pca_reference_weight),
  data.table(r, Method = "K's tau", Weight = 1 / r),
  data.table(r, Method = "Winsor", Weight = fifelse(r <= q2, 1, q2 / r)),
  data.table(r, Method = "Quad", Weight = fifelse(r <= q2, 1, q2^2 / r^2)),
  data.table(r, Method = "Ball", Weight = fifelse(r <= q2, 1, 0)),
  data.table(r, Method = "Shell", Weight = fifelse(r < q1, 0, fifelse(r <= q3, 1, 0))),
  data.table(r, Method = "LR", Weight = fifelse(
    r <= q2, 1, fifelse(r <= q3s, (q3s - r) / (q3s - q2), 0)
  ))
))

weight_order <- METHOD_ORDER[METHOD_ORDER %in% c("PCA", "K's tau", "Winsor", "Quad", "Ball", "Shell", "LR")]
weight_dt[, Method := factor(Method, levels = weight_order)]
weight_cols <- c(
  "PCA" = METHOD_COLORS[["PCA"]], "K's tau" = METHOD_COLORS[["K's tau"]],
  "Winsor" = METHOD_COLORS[["Winsor"]], "Quad" = METHOD_COLORS[["Quad"]],
  "Ball" = METHOD_COLORS[["Ball"]], "Shell" = METHOD_COLORS[["Shell"]],
  "LR" = METHOD_COLORS[["LR"]]
)
weight_types <- METHOD_LINETYPES[weight_order]
weight_labels <- c(
  "PCA" = "PCA reference",
  "K's tau" = "K's tau reference",
  "Winsor" = "Winsor",
  "Quad" = "Quad",
  "Ball" = "Ball",
  "Shell" = "Shell",
  "LR" = "LR"
)
threshold_marks <- data.table(
  r = c(q1, q2, q3, q3s),
  Weight = rep(0.075, 4),
  Label = c("Q[1]", "Q[2]", "Q[3]", "Q[3]^'*'")
)

fig02 <- ggplot(weight_dt, aes(r, Weight, colour = Method, linetype = Method)) +
  geom_line(linewidth = 0.58, lineend = "round") +
  geom_label(
    data = threshold_marks,
    aes(r, Weight, label = Label),
    inherit.aes = FALSE, parse = TRUE,
    size = 2.15, colour = "#4E5962", fill = "white",
    linewidth = 0, label.padding = unit(0.45, "mm")
  ) +
  scale_colour_manual(
    values = weight_cols, breaks = weight_order,
    labels = unname(weight_labels[weight_order])
  ) +
  scale_linetype_manual(
    values = weight_types, breaks = weight_order,
    labels = unname(weight_labels[weight_order])
  ) +
  scale_x_continuous(
    breaks = 0:5,
    limits = c(0, 5),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2), expand = c(0, 0)) +
  labs(x = "r (pairwise distance)", y = expression(xi(r)),
       colour = NULL, linetype = NULL) +
  theme_journal(6.5) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.width = unit(5.0, "mm"),
    legend.key.height = unit(2.6, "mm"),
    legend.text = element_text(size = 5.8, colour = "#30363B"),
    legend.spacing.x = unit(0.8, "mm"),
    legend.spacing.y = unit(0.2, "mm"),
    legend.margin = margin(t = 2, r = 0, b = 0, l = 0),
    plot.margin = margin(4, 5, 2, 6)
  ) +
  coord_cartesian(ylim = c(0, 2.10), clip = "on") +
  guides(
    colour = guide_legend(
      nrow = 2, byrow = TRUE,
      override.aes = list(linewidth = 0.52)
    ),
    linetype = guide_legend(nrow = 2, byrow = TRUE)
  )

weight_definition <- c(
  "PCA" = "xi(r) = 1/sqrt(2) (exact sample-covariance reference under the ordered-pair normalization)",
  "K's tau" = "xi(r) = 1/r",
  "Winsor" = "xi(r) = 1 for r <= Q2; Q2/r otherwise",
  "Quad" = "xi(r) = 1 for r <= Q2; (Q2/r)^2 otherwise",
  "Ball" = "xi(r) = 1 for r <= Q2; 0 otherwise",
  "Shell" = "xi(r) = 1 for Q1 <= r <= Q3; 0 otherwise",
  "LR" = "xi(r) = 1 for r <= Q2; linear decay to 0 at Q3*"
)
weight_dt[, `:=`(
  Definition = unname(weight_definition[as.character(Method)]),
  Threshold_values_schematic = sprintf("Q1=%.2f; Q2=%.2f; Q3=%.2f; Q3*=%.2f", q1, q2, q3, q3s),
  Visualization_range = sprintf("r from %.2f to %.2f; the right boundary is a plotting limit, not a threshold", min(r), max(r)),
  Threshold_display = "Threshold symbols label schematic transition locations; no full-height threshold guide lines are shown",
  Source_note = "Illustrative distance axis based on the original Figure 2 composition; thresholds are schematic and are not fitted experimental values. Under the ordered-pair normalization, xi(r)=1/sqrt(2) exactly reproduces the usual sample covariance matrix, whereas xi(r)=1/r recovers ordinary K's tau. The latter can exceed the displayed y-range for small r."
)]

if (!figure9_only) {
  save_drawable(
    function() print(fig02),
    file.path(out_dir, "Figure02_weight_functions"),
    width_mm = 140, height_mm = 86
  )
  weight_source <- publication_source_data(weight_dt, "r")
  # Preserve every plotted curve point while avoiding repetition of the same
  # long metadata strings on all 34,937 rows. Definitions are recorded once per
  # method and global schematic notes once for the file.
  weight_source[, Definition := fifelse(seq_len(.N) == 1L, Definition, NA_character_),
                by = Method]
  global_note_cols <- c(
    "Threshold_values_schematic", "Visualization_range",
    "Threshold_display", "Source_note"
  )
  weight_source[-1L, (global_note_cols) := NA_character_]
  fwrite(weight_source,
         file.path(source_data_dir, "Figure02_weight_functions_source_data.csv"),
         na = "")
}

if (figure2_only) {
  quit(save = "no", status = 0)
}

# -----------------------------------------------------------------------------
# Figure 9: method map and empirical evidence index
# Contract: theoretical identity is descriptive; + / +/- / - are explicitly
# qualitative summaries that direct readers to quantitative result figures.
# -----------------------------------------------------------------------------

method_summary_file <- file.path(source_data_dir, "Figure09_method_summary_source_data.csv")
if (!file.exists(method_summary_file)) stop("Missing Figure 9 summary data: ", method_summary_file)
method_summary <- fread(method_summary_file)
method_summary[, Proposed := Method %in% c("Winsor", "Quad", "Ball", "Shell", "LR")]

metric_cols <- c("Doublet / mean shift", "Dropout", "PC-wise", "Subspace", "Clustering")
score_fill <- c("+" = "#DCECE7", "+/-" = "#EEE9D7", "-" = "#E8EAEC")
score_label <- c("+" = "+", "+/-" = "±", "-" = "−")

draw_figure09 <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  x0 <- 0.025
  x1 <- 0.985
  widths <- c(0.115, 0.16, 0.16, rep(0.113, 5))
  widths <- widths / sum(widths) * (x1 - x0)
  edges <- c(x0, x0 + cumsum(widths))
  centers <- (head(edges, -1) + tail(edges, -1)) / 2

  group_y0 <- 0.925
  group_y1 <- 0.985
  group_spec <- list(
    c(1, 3, "Method characteristics", "#F1F1F1"),
    c(4, 5, "Robustness", "#F1F1F1"),
    c(6, 7, "Stability", "#F1F1F1"),
    c(8, 8, "Clustering", "#F1F1F1")
  )
  for (g in group_spec) {
    a <- as.integer(g[[1]]); b <- as.integer(g[[2]])
    grid.rect(
      x = unit((edges[a] + edges[b + 1]) / 2, "npc"),
      y = unit((group_y0 + group_y1) / 2, "npc"),
      width = unit(edges[b + 1] - edges[a], "npc"),
      height = unit(group_y1 - group_y0, "npc"),
      gp = gpar(fill = g[[4]], col = "white", lwd = 0.8)
    )
    grid.text(g[[3]], unit((edges[a] + edges[b + 1]) / 2, "npc"),
              unit((group_y0 + group_y1) / 2, "npc"),
              gp = gpar(fontfamily = "Arial", fontsize = 8.2, fontface = "bold", col = "#26323A"))
  }

  headers <- c("Method", "Core estimator", "Mechanism", metric_cols)
  header_y0 <- 0.865
  header_y1 <- group_y0
  for (j in seq_along(headers)) {
    grid.rect(x = unit(centers[j], "npc"), y = unit((header_y0 + header_y1) / 2, "npc"),
              width = unit(widths[j], "npc"), height = unit(header_y1 - header_y0, "npc"),
              gp = gpar(fill = "#FAFAFA", col = "white", lwd = 0.8))
    header_label <- if (headers[j] == "Doublet / mean shift") "Doublet /\nmean shift" else headers[j]
    grid.text(header_label, unit(centers[j], "npc"), unit((header_y0 + header_y1) / 2, "npc"),
              gp = gpar(fontfamily = "Arial", fontsize = 7.1, fontface = "bold", col = "#30363B",
                        lineheight = 0.92))
  }

  body_y0 <- 0.145
  body_y1 <- header_y0
  row_h <- (body_y1 - body_y0) / nrow(method_summary)
  for (i in seq_len(nrow(method_summary))) {
    yt <- body_y1 - (i - 1) * row_h
    yb <- yt - row_h
    yc <- (yt + yb) / 2
    char_fill <- "#FFFFFF"
    vals <- c(method_summary$Method[i], method_summary$Estimator[i], method_summary$Mechanism[i])
    if (method_summary$Proposed[i]) vals[1] <- paste0("* ", vals[1])
    for (j in 1:3) {
      grid.rect(unit(centers[j], "npc"), unit(yc, "npc"),
                width = unit(widths[j], "npc"), height = unit(row_h, "npc"),
                gp = gpar(fill = char_fill, col = "#D8DDE1", lwd = 0.45))
      grid.text(vals[j], unit(centers[j], "npc"), unit(yc, "npc"),
                gp = gpar(fontfamily = "Arial", fontsize = if (j == 1) 7.4 else 6.9,
                          fontface = if (j == 1) "bold" else "plain",
                          col = "#252A2E", lineheight = 0.92))
    }
    for (k in seq_along(metric_cols)) {
      j <- k + 3
      raw_score <- method_summary[[metric_cols[k]]][i]
      grid.rect(unit(centers[j], "npc"), unit(yc, "npc"),
                width = unit(widths[j], "npc"), height = unit(row_h, "npc"),
                gp = gpar(fill = score_fill[[raw_score]], col = "#D8DDE1", lwd = 0.45))
      grid.text(score_label[[raw_score]], unit(centers[j], "npc"), unit(yc, "npc"),
                gp = gpar(fontfamily = "Arial", fontsize = 9.4, fontface = "bold", col = "#252A2E"))
    }
  }

  evidence <- c(
    "Evidence", "", "",
    "Fig. 3;\nSupp. S1–S3",
    "Fig. 3;\nSupp. Table S3",
    "Supp. S4–S5",
    "Figs. 4–5, 7",
    "Figs. 3, 6;\nSupp. S1–S3"
  )
  ev_y0 <- 0.075
  ev_y1 <- body_y0
  for (j in seq_along(evidence)) {
    grid.rect(unit(centers[j], "npc"), unit((ev_y0 + ev_y1) / 2, "npc"),
              width = unit(widths[j], "npc"), height = unit(ev_y1 - ev_y0, "npc"),
              gp = gpar(fill = "#F5F7F9", col = "#D8DDE1", lwd = 0.45))
    grid.text(evidence[j], unit(centers[j], "npc"), unit((ev_y0 + ev_y1) / 2, "npc"),
              gp = gpar(fontfamily = "Arial", fontsize = if (j == 1) 6.4 else 6.0,
                        fontface = if (j == 1) "bold" else "plain", col = "#4B555D"))
  }
  grid.text("* Proposed methods     + relatively strong     ± mixed or condition-dependent     − relatively weak",
            unit(x0, "npc"), unit(0.038, "npc"), just = "left",
            gp = gpar(fontfamily = "Arial", fontsize = 6.1, col = "#4B555D"))
}

save_drawable(
  draw_figure09,
  file.path(out_dir, "Figure09_method_comparison"),
  width_mm = 183, height_mm = 110
)

if (figure9_only) {
  cat("Final Figure 9 outputs written to:", out_dir, "\n")
} else {
  cat("Final Figure 2 and Figure 9 outputs written to:", out_dir, "\n")
}
