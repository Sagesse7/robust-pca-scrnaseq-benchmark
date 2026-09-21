align_bottom_axis_space <- function(plots) {
  metric_file <- tempfile(fileext = ".svg")
  svglite::svglite(metric_file, width = 178 / 25.4, height = 150 / 25.4)
  on.exit({
    grDevices::dev.off()
    unlink(metric_file)
  })

  axis_height_pt <- function(p) {
    g <- ggplot2::ggplotGrob(p)
    rows <- unique(g$layout$t[grepl("^axis-b", g$layout$name)])
    max(grid::convertHeight(g$heights[rows], "pt", valueOnly = TRUE))
  }
  heights <- vapply(plots, axis_height_pt, numeric(1))
  # Retain the previous angled version's plot geometry for a controlled comparison.
  reference_heights <- vapply(plots, function(p) {
    axis_height_pt(p + ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 45, hjust = 1, vjust = 1
    )))
  }, numeric(1))
  target <- max(heights, reference_heights)
  # Compensate for differing numeric label lengths without changing their anchors.
  aligned <- Map(function(p, height) {
    p + ggplot2::theme(axis.text.x = ggplot2::element_text(
      margin = ggplot2::margin(t = 2, b = target - height, unit = "pt")
    ))
  }, plots, heights)
  checked <- vapply(aligned, axis_height_pt, numeric(1))
  stopifnot(diff(range(checked)) < 0.05)
  message("Verified equal bottom-axis space: ", round(checked[1], 2), " pt")
  aligned
}
