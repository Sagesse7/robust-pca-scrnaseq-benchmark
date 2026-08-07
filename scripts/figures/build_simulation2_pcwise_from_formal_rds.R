suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(rsvg)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: Rscript build_simulation2_pcwise_from_formal_rds.R ",
    "<formal_result_root> <authoritative_simulation2_csv> <latex_journal_root>"
  )
}

formal_root <- normalizePath(args[[1]], mustWork = TRUE)
authoritative_csv <- normalizePath(args[[2]], mustWork = TRUE)
latex_root <- normalizePath(args[[3]], mustWork = TRUE)

script_file <- sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)[1])
script_dir <- dirname(normalizePath(script_file))
release_root <- dirname(dirname(script_dir))
source(file.path(release_root, "scripts", "utils", "figure_style.R"))

source_data_dir <- file.path(latex_root, "Source_Data")
figure_dir <- file.path(latex_root, "Appendix")
dir.create(source_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

subspace_similarity <- function(A, B, k) {
  A <- qr.Q(qr(as.matrix(A)[, seq_len(min(k, ncol(A))), drop = FALSE]))
  B <- qr.Q(qr(as.matrix(B)[, seq_len(min(k, ncol(B))), drop = FALSE]))
  singular_values <- svd(crossprod(A, B), nu = 0, nv = 0)$d
  mean(pmin(pmax(singular_values, 0), 1))
}

pcwise_absolute_cosine <- function(A, B, k = 20L) {
  k <- min(k, ncol(A), ncol(B))
  A <- as.matrix(A)[, seq_len(k), drop = FALSE]
  B <- as.matrix(B)[, seq_len(k), drop = FALSE]
  denominator <- sqrt(colSums(A^2) * colSums(B^2))
  values <- abs(colSums(A * B) / denominator)
  pmin(pmax(values, 0), 1)
}

resolve_formal_path <- function(relative_path) {
  file.path(formal_root, as.character(relative_path))
}

sim2 <- fread(authoritative_csv)
required_columns <- c(
  "BaseID", "Replicate", "Method", "Alpha", "ParamID", "Condition",
  "dropout_mid", "doublet_frac", "cell_frac", "meanlog", "GeneHash",
  "CleanPath", "RDS_Path", "SimilarityPC10", "SimilarityPC20"
)
missing_columns <- setdiff(required_columns, names(sim2))
if (length(missing_columns)) {
  stop("Authoritative CSV is missing columns: ", paste(missing_columns, collapse = ", "))
}

sim2 <- sim2[is.na(Alpha) | abs(Alpha - 0.05) < 1e-12]
if (nrow(sim2) != 2000L) {
  stop("Expected 2,000 alpha=0.05 Simulation 2 rows; found ", nrow(sim2))
}

sim2[, `:=`(
  CleanResolved = resolve_formal_path(CleanPath),
  PerturbedResolved = resolve_formal_path(RDS_Path)
)]
missing_files <- unique(c(
  sim2[!file.exists(CleanResolved), CleanResolved],
  sim2[!file.exists(PerturbedResolved), PerturbedResolved]
))
if (length(missing_files)) {
  stop("Missing formal RDS files (first shown): ", missing_files[[1]])
}

clean_cache <- new.env(parent = emptyenv())
get_clean <- function(path) {
  if (!exists(path, envir = clean_cache, inherits = FALSE)) {
    assign(path, readRDS(path), envir = clean_cache)
  }
  get(path, envir = clean_cache, inherits = FALSE)
}

pcwise_rows <- vector("list", nrow(sim2))
validation_rows <- vector("list", nrow(sim2))

for (i in seq_len(nrow(sim2))) {
  row <- sim2[i]
  clean_obj <- get_clean(row$CleanResolved)
  perturbed_obj <- readRDS(row$PerturbedResolved)

  clean_method_name <- publication_method(clean_obj$Method)
  perturbed_method_name <- publication_method(perturbed_obj$Method)
  expected_method <- publication_method(row$Method)

  metadata_ok <- identical(as.integer(clean_obj$BaseID), as.integer(row$BaseID)) &&
    identical(as.integer(perturbed_obj$BaseID), as.integer(row$BaseID)) &&
    identical(as.integer(clean_obj$Replicate), as.integer(row$Replicate)) &&
    identical(as.integer(perturbed_obj$Replicate), as.integer(row$Replicate)) &&
    identical(as.character(clean_method_name), as.character(expected_method)) &&
    identical(as.character(perturbed_method_name), as.character(expected_method)) &&
    identical(as.integer(perturbed_obj$ParamID), as.integer(row$ParamID)) &&
    identical(clean_obj$gene_universe, perturbed_obj$gene_universe)

  if (!metadata_ok) {
    stop("RDS metadata mismatch at authoritative row ", i)
  }

  clean_loadings <- clean_obj$loadings_full
  perturbed_loadings <- perturbed_obj$loadings_full
  if (nrow(clean_loadings) != nrow(perturbed_loadings) ||
      ncol(clean_loadings) < 20L || ncol(perturbed_loadings) < 20L) {
    stop("Incompatible loading matrices at authoritative row ", i)
  }

  pc_cosine <- pcwise_absolute_cosine(clean_loadings, perturbed_loadings, 20L)
  recomputed_pc10 <- subspace_similarity(clean_loadings, perturbed_loadings, 10L)
  recomputed_pc20 <- subspace_similarity(clean_loadings, perturbed_loadings, 20L)

  pcwise_rows[[i]] <- data.table(
    BaseID = row$BaseID,
    Replicate = row$Replicate,
    ParamID = row$ParamID,
    Condition = row$Condition,
    dropout_mid = row$dropout_mid,
    doublet_frac = row$doublet_frac,
    cell_frac = row$cell_frac,
    meanlog = row$meanlog,
    Method = expected_method,
    PC = seq_len(20L),
    AbsoluteCosine = pc_cosine
  )

  validation_rows[[i]] <- data.table(
    BaseID = row$BaseID,
    Replicate = row$Replicate,
    ParamID = row$ParamID,
    Method = expected_method,
    AuthoritativeSimilarityPC10 = row$SimilarityPC10,
    RecomputedSimilarityPC10 = recomputed_pc10,
    DifferencePC10 = recomputed_pc10 - row$SimilarityPC10,
    AuthoritativeSimilarityPC20 = row$SimilarityPC20,
    RecomputedSimilarityPC20 = recomputed_pc20,
    DifferencePC20 = recomputed_pc20 - row$SimilarityPC20,
    MetadataMatch = metadata_ok
  )

  if (i %% 100L == 0L) {
    message("Processed ", i, " / ", nrow(sim2), " paired formal results")
  }
}

pcwise_long <- rbindlist(pcwise_rows)
validation <- rbindlist(validation_rows)
max_difference <- max(abs(c(validation$DifferencePC10, validation$DifferencePC20)))
if (!is.finite(max_difference) || max_difference > 1e-12) {
  stop("Recomputed subspace similarity differs from the authoritative CSV; max |difference| = ", max_difference)
}

replicate_level <- pcwise_long[, .(
  PCwisePC10 = mean(AbsoluteCosine[PC <= 10L]),
  PCwisePC20 = mean(AbsoluteCosine[PC <= 20L])
), by = .(
  BaseID, Replicate, ParamID, Condition, dropout_mid, doublet_frac,
  cell_frac, meanlog, Method
)]

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

sweeps <- make_sweeps(replicate_level)
figure_source <- rbindlist(lapply(names(sweeps), function(key) {
  spec <- sweeps[[key]]
  z <- copy(spec$data)
  z[, `:=`(
    Noise = spec$xlab,
    Level = factor(level, levels = spec$levels, labels = spec$labels)
  )]
  z[, .(
    Mean = mean(PCwisePC10),
    SD = sd(PCwisePC10),
    Median = median(PCwisePC10),
    Q1 = quantile(PCwisePC10, 0.25, names = FALSE),
    Q3 = quantile(PCwisePC10, 0.75, names = FALSE),
    N = .N
  ), by = .(Noise, Level, Method)]
}))

publication_long <- copy(pcwise_long)
setorder(publication_long, ParamID, Replicate, Method, PC)
fwrite(
  publication_long,
  file.path(source_data_dir, "FigureS04_simulation2_pcwise_pc1_pc20_source_data.csv")
)

publication_replicate <- copy(replicate_level)
setorder(publication_replicate, ParamID, Replicate, Method)
fwrite(
  publication_replicate,
  file.path(source_data_dir, "FigureS04_simulation2_pcwise_pc10_replicate_source_data.csv")
)

publication_summary <- copy(figure_source)
setnames(publication_summary, "Noise", "PerturbationCondition")
publication_summary[, MethodOrder := match(Method, PUBLICATION_METHOD_ORDER)]
setorder(publication_summary, PerturbationCondition, Level, MethodOrder)
publication_summary[, MethodOrder := NULL]
fwrite(
  publication_summary,
  file.path(source_data_dir, "FigureS04_simulation2_pcwise_pc10_source_data.csv")
)

fwrite(
  validation,
  file.path(source_data_dir, "FigureS04_simulation2_pcwise_validation.csv")
)

distribution_panel <- function(spec, key) {
  z <- copy(spec$data[is.finite(PCwisePC10)])
  z[, Method := factor(as.character(clean_method(Method)), levels = METHOD_ORDER)]
  z[, Level := factor(level, levels = spec$levels, labels = spec$labels)]
  reference <- z[0]
  if (identical(key, "dropout")) {
    reference <- z[level == "0", .(PCwisePC10 = mean(PCwisePC10)), by = Method]
    reference[, Level := factor("0", levels = spec$labels)]
    z <- z[level != "0"]
  }

  axis_labels <- spec$labels
  axis_title <- spec$xlab
  if (key %in% c("doublet", "shift_fraction")) {
    axis_labels <- sub("%$", "", spec$labels)
    axis_title <- paste0(spec$xlab, " (%)")
  }

  jitter_seed <- c(
    dropout = 451L, doublet = 452L,
    shift_fraction = 453L, shift_factor = 454L
  )[[key]]

  p <- ggplot(z, aes(Level, PCwisePC10)) +
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
        "rect", xmin = reference_position - 0.5, xmax = reference_position + 0.5,
        ymin = -Inf, ymax = Inf, fill = "#E2E2E2", alpha = 0.55
      ) +
      geom_point(
        data = reference, aes(Level, PCwisePC10), inherit.aes = FALSE,
        shape = 21, fill = "white", colour = "#555555", size = 1.05, stroke = 0.3
      )
  }
  p
}

panels <- Map(distribution_panel, sweeps, names(sweeps))
figure_s4 <- wrap_plots(panels, ncol = 2) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

save_pub_r(
  figure_s4,
  file.path(figure_dir, "FigureS04_simulation2_pcwise_pc10"),
  183, 156
)

rank_panel <- function(spec, key) {
  z <- copy(spec$data[is.finite(PCwisePC10)])
  z[, Method := factor(as.character(clean_method(Method)), levels = METHOD_ORDER)]
  z[, Level := factor(level, levels = spec$levels, labels = spec$labels)]

  ranks <- z[, .(
    Median = median(PCwisePC10),
    N = .N
  ), by = .(Method, Level)]
  ranks[, Rank := frank(-Median, ties.method = "average"), by = Level]

  if (identical(key, "dropout")) {
    ranks[as.character(Level) == "0", Rank := NA_real_]
  }

  ranks[, `:=`(
    Label = ifelse(is.na(Rank), "NR", sprintf("%.0f", Rank)),
    TextColour = ifelse(!is.na(Rank) & Rank >= 6, "white", "#151515"),
    Noise = spec$xlab
  )]

  source_ranks <- copy(ranks)
  source_ranks[, `:=`(
    LevelOrder = as.integer(Level),
    MethodOrder = as.integer(Method),
    Method = publication_method(Method),
    Level = as.character(Level)
  )]
  setorder(source_ranks, LevelOrder, MethodOrder)

  plot <- ggplot(ranks, aes(Level, Method, fill = Rank)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(aes(label = Label, colour = TextColour), size = 1.8) +
    scale_colour_identity() +
    scale_y_discrete(
      limits = rev(METHOD_ORDER),
      labels = METHOD_LABELS
    ) +
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

  list(plot = plot, source = source_ranks)
}

rank_results <- Map(rank_panel, sweeps, names(sweeps))
rank_source <- rbindlist(lapply(rank_results, `[[`, "source"))
rank_source[, NoiseOrder := match(
  Noise,
  c("Dropout midpoint", "Doublet fraction", "Shifted-cell fraction", "Shift magnitude")
)]
setorder(rank_source, NoiseOrder, LevelOrder, MethodOrder)
rank_source[, c("NoiseOrder", "LevelOrder", "MethodOrder", "TextColour") := NULL]
setnames(rank_source, "Noise", "PerturbationCondition")
fwrite(
  rank_source,
  file.path(source_data_dir, "FigureS05_simulation2_pcwise_rank_source_data.csv")
)

figure_s5 <- wrap_plots(lapply(rank_results, `[[`, "plot"), ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 8, face = "bold"))

save_pub_r(
  figure_s5,
  file.path(figure_dir, "FigureS05_simulation2_pcwise_rank"),
  183, 112
)

cat("Simulation 2 PC-wise figures generated.\n")
cat("Paired formal rows:", nrow(replicate_level), "\n")
cat("Per-PC rows:", nrow(pcwise_long), "\n")
cat("Unique clean RDS files:", length(ls(clean_cache)), "\n")
cat("Maximum authoritative subspace validation difference:", format(max_difference, digits = 16), "\n")
