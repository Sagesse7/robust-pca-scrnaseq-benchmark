parse_cli <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) next
    p <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    out[[p[1]]] <- paste(p[-1], collapse = "=")
  }
  out
}

cli_csv <- function(x, default) {
  if (is.null(x) || !nzchar(x)) return(default)
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

cli_numeric_csv <- function(x, default) as.numeric(cli_csv(x, as.character(default)))

to_int_labels <- function(x) {
  if (is.null(x)) stop("Labels are NULL.")
  if (is.numeric(x) || is.integer(x)) return(as.integer(x))
  as.integer(factor(as.character(x)))
}

normalize_log_scale_counts <- function(C) {
  C <- as.matrix(C)
  lib <- colSums(C)
  lib[lib == 0] <- 1
  X <- t(log1p(t(1e6 * t(C) / lib)))
  keep <- apply(X, 2, function(z) isTRUE(is.finite(var(z))) && var(z) > 0)
  X <- X[, keep, drop = FALSE]
  X <- X[, apply(X, 2, function(z) all(is.finite(z))), drop = FALSE]
  X <- scale(X, center = TRUE, scale = TRUE)
  X[, apply(X, 2, function(z) all(is.finite(z))), drop = FALSE]
}

add_gene_subset_shift <- function(C, cell_frac, gene_frac, fold_change, sdlog) {
  if (is.na(cell_frac) || is.na(gene_frac) || is.na(fold_change) || is.na(sdlog) ||
      cell_frac <= 0 || gene_frac <= 0) {
    return(list(C = C, cells = integer()))
  }
  cells <- sample.int(ncol(C), max(1L, round(cell_frac * ncol(C))))
  ng <- max(1L, round(gene_frac * nrow(C)))
  for (j in cells) {
    genes <- sample.int(nrow(C), ng)
    C[genes, j] <- C[genes, j] * rlnorm(1L, log(fold_change), sdlog)
  }
  list(C = C, cells = cells)
}

add_doublets <- function(C, fraction) {
  if (is.na(fraction) || fraction <= 0) return(list(C = C, cells = integer()))
  source_counts <- C
  cells <- sample.int(ncol(C), max(1L, round(fraction * ncol(C))))
  for (j in cells) {
    donors <- sample(setdiff(seq_len(ncol(C)), j), 2L)
    C[, j] <- source_counts[, donors[1]] + source_counts[, donors[2]]
  }
  list(
    C = C,
    cells = cells,
    donor_design = "two_unmodified_background_singlets"
  )
}

apply_splatter_dropout <- function(sce_base, midpoint, shape = NULL) {
  if (is.na(midpoint)) stop("A finite dropout midpoint is required.")
  params <- splatter::newSplatParams()
  params <- splatter::setParam(params, "nGenes", nrow(sce_base))
  params <- splatter::setParam(params, "batchCells", ncol(sce_base))
  params <- splatter::setParam(params, "dropout.type", "experiment")
  params <- splatter::setParam(params, "dropout.mid", midpoint)
  if (!is.null(shape)) params <- splatter::setParam(params, "dropout.shape", shape)
  getFromNamespace("splatSimDropout", "splatter")(sce_base, params)
}

simulate_clean_splatter_base <- function(p, seed, de_fac_loc = 0.1) {
  set.seed(seed)
  gp <- seq_len(p$CellType)
  gp <- gp / sum(gp)
  sp <- splatter::newSplatParams(batchCells = p$Cell, group.prob = gp)
  sp <- splatter::setParam(sp, "nGenes", p$Gene)
  sp <- splatter::setParam(sp, "de.prob", p$de_prob)
  sp <- splatter::setParam(sp, "de.facLoc", de_fac_loc)
  splatter::splatSimulate(
    sp, method = "groups", dropout.type = "none", verbose = FALSE
  )
}

isolated_condition_name <- function(p) {
  active <- c(
    dropout = !is.na(p$dropout_mid),
    doublet = !is.na(p$doublet_frac) && p$doublet_frac > 0,
    shift = !is.na(p$cell_frac) && p$cell_frac > 0
  )
  if (sum(active) > 1L) stop("Isolated-noise design requires at most one active perturbation.")
  if (!any(active)) "clean" else names(active)[active]
}

condition_seed_offset <- function(condition) {
  switch(
    condition,
    clean = 0L,
    shift = 11L,
    doublet = 22L,
    dropout = 33L,
    stop("Unsupported condition: ", condition)
  )
}

perturbation_seed <- function(clean_seed, condition) {
  if (condition == "clean") return(NA_integer_)
  as.integer(clean_seed + condition_seed_offset(condition))
}

background_dropout_seed <- function(clean_seed) as.integer(clean_seed + 501L)

make_isolated_condition <- function(sce_base, p, seed) {
  condition <- isolated_condition_name(p)
  condition_seed <- perturbation_seed(seed, condition)
  base_counts <- as.matrix(SummarizedExperiment::assay(sce_base, "counts"))
  C <- base_counts
  contaminated <- rep(FALSE, ncol(C))
  doublet_donor_design <- NA_character_

  if (condition == "dropout") {
    set.seed(condition_seed)
    dropped <- apply_splatter_dropout(sce_base, p$dropout_mid)
    C <- as.matrix(SummarizedExperiment::assay(dropped, "counts"))
  } else if (condition == "doublet") {
    set.seed(condition_seed)
    z <- add_doublets(C, p$doublet_frac)
    C <- z$C
    contaminated[z$cells] <- TRUE
    doublet_donor_design <- z$donor_design
  } else if (condition == "shift") {
    set.seed(condition_seed)
    z <- add_gene_subset_shift(C, p$cell_frac, p$gene_frac, p$meanlog, p$sdlog)
    C <- z$C
    contaminated[z$cells] <- TRUE
  }

  list(
    condition = condition,
    base_counts = base_counts,
    counts = C,
    labels = to_int_labels(SummarizedExperiment::colData(sce_base)$Group),
    contaminated = contaminated,
    doublet_donor_design = doublet_donor_design,
    perturbation_seed = condition_seed
  )
}

make_background_dropout_condition <- function(sce_base, p, seed, background_mid = 0) {
  condition <- isolated_condition_name(p)
  condition_seed <- perturbation_seed(seed, condition)
  background_seed <- background_dropout_seed(seed)
  original_base_counts <- as.matrix(SummarizedExperiment::assay(sce_base, "counts"))

  set.seed(background_seed)
  background_sce <- apply_splatter_dropout(sce_base, background_mid)
  background_counts <- as.matrix(SummarizedExperiment::assay(background_sce, "counts"))

  C <- background_counts
  contaminated <- rep(FALSE, ncol(C))
  dropout_seed <- NA_integer_
  doublet_donor_design <- NA_character_

  if (condition == "dropout") {
    # The dropout sweep directly sets the total dropout level. It is not stacked
    # on top of the common background dropout.
    dropout_seed <- background_seed
    set.seed(dropout_seed)
    dropped <- apply_splatter_dropout(sce_base, p$dropout_mid)
    C <- as.matrix(SummarizedExperiment::assay(dropped, "counts"))
  } else if (condition == "doublet") {
    set.seed(condition_seed)
    z <- add_doublets(C, p$doublet_frac)
    C <- z$C
    contaminated[z$cells] <- TRUE
    doublet_donor_design <- z$donor_design
  } else if (condition == "shift") {
    set.seed(condition_seed)
    z <- add_gene_subset_shift(C, p$cell_frac, p$gene_frac, p$meanlog, p$sdlog)
    C <- z$C
    contaminated[z$cells] <- TRUE
  }

  list(
    condition = condition,
    original_base_counts = original_base_counts,
    background_counts = background_counts,
    counts = C,
    labels = to_int_labels(SummarizedExperiment::colData(sce_base)$Group),
    contaminated = contaminated,
    doublet_donor_design = doublet_donor_design,
    background_mid = background_mid,
    background_seed = background_seed,
    dropout_seed = dropout_seed,
    perturbation_seed = condition_seed
  )
}

simulation_condition_diagnostics <- function(
    base_counts, condition_counts, condition, contaminated, analysis_genes,
    labels = NULL
) {
  base_counts <- as.matrix(base_counts)
  condition_counts <- as.matrix(condition_counts)
  originally_nonzero <- base_counts != 0
  newly_zero <- originally_nonzero & condition_counts == 0
  out <- list(
    NCells = ncol(condition_counts),
    InputGenes = nrow(condition_counts),
    AnalysisGenes = analysis_genes,
    DOverN = analysis_genes / ncol(condition_counts),
    BaseZeroFraction = mean(base_counts == 0),
    ConditionZeroFraction = mean(condition_counts == 0),
    AddedZeroFractionAllEntries = mean(newly_zero),
    AddedZeroFractionAmongNonzero = if (any(originally_nonzero)) {
      sum(newly_zero) / sum(originally_nonzero)
    } else NA_real_,
    ChangedEntryFraction = mean(base_counts != condition_counts),
    ContaminatedCells = sum(contaminated),
    ContaminatedCellFraction = mean(contaminated),
    ConditionDiagnostic = condition
  )
  if (!is.null(labels)) {
    group_counts <- table(to_int_labels(labels))
    out$GroupCountMin <- min(group_counts)
    out$GroupCountMax <- max(group_counts)
    out$GroupCountVector <- paste(
      paste0(names(group_counts), ":", as.integer(group_counts)),
      collapse = ";"
    )
  }
  out
}

set_single_thread_blas <- function() {
  Sys.setenv(
    OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1"
  )
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
  }
}

nonempty_output_guard <- function(path) {
  if (dir.exists(path) && length(list.files(path, all.files = TRUE, no.. = TRUE))) {
    stop("Output directory is not empty; refusing to overwrite: ", path)
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}
