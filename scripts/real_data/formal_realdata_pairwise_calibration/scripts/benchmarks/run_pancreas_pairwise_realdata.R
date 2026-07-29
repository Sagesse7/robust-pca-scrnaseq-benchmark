# =========================================================
# Real-data benchmark (parallel; Simulation1/2 style)
# FULL data as reference + repeated random subset(2000 cells)
# Compare FULL vs SUBSET clustering by ARI / NMI
# Also compare FULL / SUBSET against true labels
# Also compare FULL vs SUBSET PC similarity in gene space
# Select top variable genes before scaling
# Repeat clustering (k-means / Louvain / GMM) multiple times
# =========================================================

suppressPackageStartupMessages({
  library(mclust)      # ARI, Mclust(GMM)
  library(aricode)     # NMI
  library(rrcov)       # PcaHubert, PcaGrid
  library(scran)       # buildSNNGraph
  library(igraph)      # Louvain
  library(rpca)        # PCP
  
  library(future)
  library(future.apply)
})

# -----------------------------
# User settings
# -----------------------------
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  getwd()
}

SCRIPT_DIR <- get_script_dir()
PROJECT_DIR <- dirname(dirname(SCRIPT_DIR))

parse_cli <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--") || !grepl("=", arg, fixed = TRUE)) next
    parts <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    out[[parts[1]]] <- paste(parts[-1], collapse = "=")
  }
  out
}

parse_methods <- function(x) unique(trimws(strsplit(x, ",", fixed = TRUE)[[1]]))
cli <- parse_cli(commandArgs(trailingOnly = TRUE))
DATA_ROOT <- if (is.null(cli[["data-root"]])) file.path(PROJECT_DIR, "data") else cli[["data-root"]]
CSV_PATH   <- file.path(DATA_ROOT, "pancreas_GSM2230759_human3_umifm_counts.csv")
LABEL_COL  <- "assigned_cluster"
run_id <- paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"))
OUT_DIR    <- if (is.null(cli[["output-dir"]])) file.path(PROJECT_DIR, "results", "pancreas", run_id) else cli[["output-dir"]]
OUT_CSV    <- file.path(OUT_DIR, "realdata_full_vs_subset_results_with_pc_similarity.csv")

SUBSET_N   <- if (is.null(cli[["subset-n"]])) 2000L else as.integer(cli[["subset-n"]])
REPS       <- if (is.null(cli$repeats)) 20L else as.integer(cli$repeats)
CLUSTER_REPS <- if (is.null(cli[["cluster-repeats"]])) 10L else as.integer(cli[["cluster-repeats"]])
BASE_SEED  <- if (is.null(cli[["base-seed"]])) 12345L else as.integer(cli[["base-seed"]])
MAX_PC     <- 20L

RUN_GMM     <- TRUE
RUN_LOUVAIN <- TRUE
RUN_KMEANS  <- TRUE

LOUVAIN_K  <- 15L
KMAX_MODEL <- 15L

WORKERS <- if (is.null(cli$workers)) 5L else as.integer(cli$workers)

TOP_HVG <- if (is.null(cli$hvg)) 2000L else as.integer(cli$hvg)
PAIRWISE_ALPHA <- if (is.null(cli$alpha)) 0.05 else as.numeric(cli$alpha)

METHOD_NAMES <- parse_methods(if (is.null(cli$methods)) {
  "PCA,Grid,Hubert,PCP,Winsor,Quad,Ball,Shell,LR,Tau"
} else cli$methods)
SUPPORTED_METHODS <- c("PCA", "Grid", "Hubert", "PCP", "Winsor", "Quad", "Ball", "Shell", "LR", "Tau")
if (!all(METHOD_NAMES %in% SUPPORTED_METHODS)) stop("Unsupported methods: ", paste(setdiff(METHOD_NAMES, SUPPORTED_METHODS), collapse = ","))
if (!is.finite(PAIRWISE_ALPHA) || PAIRWISE_ALPHA <= 0 || PAIRWISE_ALPHA >= 1) stop("--alpha must lie strictly between 0 and 1.")

# -----------------------------
# HPC safety: prevent nested threading
# -----------------------------
suppressPackageStartupMessages({
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
  }
  Sys.setenv(
    "OPENBLAS_NUM_THREADS"   = "1",
    "MKL_NUM_THREADS"        = "1",
    "OMP_NUM_THREADS"        = "1",
    "VECLIB_MAXIMUM_THREADS" = "1"
  )
})

# -----------------------------
# Basic helpers
# -----------------------------
to_int <- function(v) {
  if (is.factor(v)) return(as.integer(v))
  if (is.character(v)) return(as.integer(factor(v)))
  if (is.numeric(v) || is.integer(v)) return(as.integer(v))
  as.integer(factor(as.character(v)))
}

ari_safe <- function(x, y) {
  if (is.null(x) || is.null(y)) return(NA_real_)
  if (length(x) != length(y) || length(x) == 0) return(NA_real_)
  tryCatch(mclust::adjustedRandIndex(x, y), error = function(e) NA_real_)
}

nmi_safe <- function(x, y) {
  if (is.null(x) || is.null(y)) return(NA_real_)
  if (length(x) != length(y) || length(x) == 0) return(NA_real_)
  tryCatch(aricode::NMI(x, y), error = function(e) NA_real_)
}

safe_var_keep <- function(X) {
  keep <- apply(X, 2, function(z) var(z, na.rm = TRUE) > 0)
  keep[is.na(keep)] <- FALSE
  X[, keep, drop = FALSE]
}

safe_scale <- function(X) {
  Xs <- scale(X, center = TRUE, scale = TRUE)
  Xs <- Xs[, colSums(!is.finite(Xs)) == 0, drop = FALSE]
  Xs
}

select_pc <- function(emb, p) {
  if (is.null(emb)) return(NULL)
  if (!is.matrix(emb)) emb <- as.matrix(emb)
  k <- min(p, ncol(emb))
  if (k <= 0) return(NULL)
  emb[, seq_len(k), drop = FALSE]
}

select_loading <- function(V, p) {
  if (is.null(V)) return(NULL)
  V <- as.matrix(V)
  k <- min(p, ncol(V))
  if (k <= 0) return(NULL)
  V[, seq_len(k), drop = FALSE]
}

normalize_columns <- function(V) {
  V <- as.matrix(V)
  cn <- sqrt(colSums(V^2))
  cn[cn == 0 | !is.finite(cn)] <- 1
  sweep(V, 2, cn, "/", check.margin = FALSE)
}

pcwise_abs_cosine <- function(V_full, V_sub, max_pc = 20L) {
  out <- setNames(rep(NA_real_, max_pc), paste0("pc_cos_", seq_len(max_pc)))
  if (is.null(V_full) || is.null(V_sub)) return(out)
  
  V1 <- normalize_columns(V_full)
  V2 <- normalize_columns(V_sub)
  
  k <- min(max_pc, ncol(V1), ncol(V2), nrow(V1), nrow(V2))
  if (k <= 0) return(out)
  
  sims <- colSums(V1[, seq_len(k), drop = FALSE] * V2[, seq_len(k), drop = FALSE])
  out[seq_len(k)] <- abs(sims)
  out
}

subspace_similarity <- function(V_full, V_sub, p = 20L) {
  out <- list(
    mean_cos = NA_real_,
    min_cos  = NA_real_,
    frob     = NA_real_
  )
  
  if (is.null(V_full) || is.null(V_sub)) return(out)
  
  U <- select_loading(V_full, p)
  V <- select_loading(V_sub, p)
  if (is.null(U) || is.null(V)) return(out)
  
  k <- min(ncol(U), ncol(V))
  if (k <= 0) return(out)
  
  U <- U[, seq_len(k), drop = FALSE]
  V <- V[, seq_len(k), drop = FALSE]
  
  Uq <- qr.Q(qr(U))
  Vq <- qr.Q(qr(V))
  
  sv <- tryCatch(
    svd(t(Uq) %*% Vq, nu = 0, nv = 0)$d,
    error = function(e) rep(NA_real_, k)
  )
  
  sv <- pmin(1, pmax(0, sv))
  
  out$mean_cos <- mean(sv, na.rm = TRUE)
  out$min_cos  <- min(sv, na.rm = TRUE)
  out$frob     <- sqrt(sum(sv^2, na.rm = TRUE) / k)
  out
}

append_subspace_prefix <- function(x, prefix) {
  names(x) <- paste0(prefix, names(x))
  x
}

append_prefix <- function(lst, prefix) {
  names(lst) <- paste0(prefix, names(lst))
  lst
}

cutoff_columns <- function(cutoff, prefix) {
  fields <- c("Q1", "Q2", "Q3", "Q3star", "mean_weight", "median_weight",
              "prop_weight_lt_1", "prop_weight_lt_0.5", "prop_weight_zero",
              "effective_pairs", "total_pairs")
  values <- setNames(as.list(rep(NA_real_, length(fields))), paste0(prefix, fields))
  if (!is.null(cutoff)) {
    for (field in fields) values[[paste0(prefix, field)]] <- cutoff[[field]]
  }
  values
}

# -----------------------------
# Read real CSV
# -----------------------------
read_real_csv_pancreas <- function(csv_path, label_col = "assigned_cluster") {
  dat <- read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE)
  
  if (!(label_col %in% colnames(dat))) {
    stop(sprintf("Column '%s' not found in CSV.", label_col))
  }
  
  labels <- dat[[label_col]]
  
  id_col <- NULL
  if (colnames(dat)[1] == "") {
    id_col <- dat[[1]]
  }
  
  meta_cols <- c("", "barcode", label_col)
  expr_cols <- setdiff(colnames(dat), meta_cols)
  
  is_num <- vapply(dat[, expr_cols, drop = FALSE], is.numeric, logical(1))
  expr_cols <- expr_cols[is_num]
  
  if (length(expr_cols) < 2) {
    stop("Fewer than 2 numeric expression columns remain after excluding metadata.")
  }
  
  X <- as.matrix(dat[, expr_cols, drop = FALSE])   # cells x genes
  storage.mode(X) <- "numeric"
  
  if (!is.null(id_col)) {
    rownames(X) <- make.unique(as.character(id_col))
  } else {
    rownames(X) <- sprintf("cell_%07d", seq_len(nrow(X)))
  }
  
  list(
    counts = X,
    labels = labels,
    genes  = expr_cols
  )
}

# -----------------------------
# Normalize + log + top-HVG + global scale
# -----------------------------
preprocess_real_counts <- function(C_cells_genes, top_hvg = 2000L) {
  C <- as.matrix(C_cells_genes)
  storage.mode(C) <- "numeric"
  
  keep_gene <- colSums(C, na.rm = TRUE) > 0
  C <- C[, keep_gene, drop = FALSE]
  
  libsize <- rowSums(C, na.rm = TRUE)
  libsize[libsize == 0] <- 1
  
  CPM <- C / libsize * 1e6
  LOGCPM <- log1p(CPM)
  
  LOGCPM <- safe_var_keep(LOGCPM)
  
  gene_var <- apply(LOGCPM, 2, var, na.rm = TRUE)
  gene_var[is.na(gene_var)] <- -Inf
  
  k_hvg <- min(as.integer(top_hvg), ncol(LOGCPM))
  if (k_hvg <= 0) stop("No genes available after variance filtering.")
  
  ord <- order(gene_var, decreasing = TRUE)
  keep_hvg <- ord[seq_len(k_hvg)]
  LOGCPM <- LOGCPM[, keep_hvg, drop = FALSE]
  
  X <- safe_scale(LOGCPM)
  
  list(
    X = X,
    selected_genes = colnames(X),
    n_hvg = ncol(X)
  )
}

# -----------------------------
# PCP
# -----------------------------
pcp_decompose_fast <- function(X, lambda = 1 / sqrt(max(dim(X)))) {
  res <- rpca::rpca(
    X,
    lambda     = lambda,
    term.delta = 1e-5,
    max.iter   = 500,
    trace      = FALSE
  )
  list(L = res$L, S = res$S)
}

# -----------------------------
# PcaGrid safe wrapper
# -----------------------------
pca_grid_fit_safe <- function(X, k = 20) {
  X <- as.matrix(X)
  k <- min(k, nrow(X), ncol(X))
  p <- ncol(X)
  
  fit <- rrcov::PcaGrid(
    as.data.frame(X),
    k = k, kmax = k,
    center = rep(0, p),
    scale  = rep(1, p)
  )
  
  list(
    scores   = as.matrix(fit@scores)[, seq_len(k), drop = FALSE],
    loadings = as.matrix(fit@loadings)[, seq_len(k), drop = FALSE]
  )
}

# -----------------------------
# Kendall and generalized Kendall methods
# -----------------------------
# The canonical implementations, including the numerically stable K's tau
# estimator, are shared across all three real-data scripts.
source(file.path(SCRIPT_DIR, "pairwise_calibration_core.R"))

# -----------------------------
# Methods builder
# -----------------------------
build_methods <- function(max_pc = 20) {
  list(
    PCA = function(X) {
      k <- min(max_pc, nrow(X), ncol(X))
      fit <- prcomp(X, center = FALSE, scale. = FALSE)
      list(
        scores   = as.matrix(fit$x)[, seq_len(k), drop = FALSE],
        loadings = as.matrix(fit$rotation)[, seq_len(k), drop = FALSE]
      )
    },
    Grid = function(X) {
      pca_grid_fit_safe(X, k = max_pc)
    },
    Hubert = function(X) {
      k <- min(max_pc, nrow(X), ncol(X))
      fit <- rrcov::PcaHubert(X, k = k, kmax = k, scale = FALSE)
      list(
        scores   = as.matrix(fit@scores)[, seq_len(k), drop = FALSE],
        loadings = as.matrix(fit@loadings)[, seq_len(k), drop = FALSE]
      )
    },
    PCP = function(X) {
      k <- min(max_pc, nrow(X), ncol(X))
      rp <- pcp_decompose_fast(X)
      fit <- prcomp(rp$L, center = FALSE, scale. = FALSE)
      V <- as.matrix(fit$rotation)[, seq_len(k), drop = FALSE]
      list(
        scores   = X %*% V,
        loadings = V
      )
    },
    Winsor = function(X) {
      k <- min(max_pc, nrow(X), ncol(X))
      eg <- function_winsor_vectorized(X)
      V  <- as.matrix(eg$vectors)[, seq_len(k), drop = FALSE]
      list(scores = X %*% V, loadings = V, cutoff = eg$cutoff)
    },
    Quad = function(X) {
      k <- min(max_pc, nrow(X), ncol(X))
      eg <- function_quad_vectorized(X)
      V  <- as.matrix(eg$vectors)[, seq_len(k), drop = FALSE]
      list(scores = X %*% V, loadings = V, cutoff = eg$cutoff)
    },
    Ball = function(X) {
      k <- min(max_pc, nrow(X), ncol(X))
      eg <- function_ball_vectorized(X)
      V  <- as.matrix(eg$vectors)[, seq_len(k), drop = FALSE]
      list(scores = X %*% V, loadings = V, cutoff = eg$cutoff)
    },
    Shell = function(X) {
      k <- min(max_pc, nrow(X), ncol(X))
      eg <- function_shell_vectorized(X)
      V  <- as.matrix(eg$vectors)[, seq_len(k), drop = FALSE]
      list(scores = X %*% V, loadings = V, cutoff = eg$cutoff)
    },
    LR = function(X) {
      k <- min(max_pc, nrow(X), ncol(X))
      eg <- function_lr_vectorized(X)
      V  <- as.matrix(eg$vectors)[, seq_len(k), drop = FALSE]
      list(scores = X %*% V, loadings = V, cutoff = eg$cutoff)
    },
    Tau = function(X) {
      k <- min(max_pc, nrow(X), ncol(X))
      eg <- function_tau_vectorized(X)
      V  <- as.matrix(eg$vectors)[, seq_len(k), drop = FALSE]
      list(scores = X %*% V, loadings = V)
    }
  )
}

safe_method_run <- function(method_fun, X, method_name = "unknown") {
  tryCatch({
    res <- method_fun(X)
    
    scores <- as.matrix(res$scores)
    loadings <- as.matrix(res$loadings)
    
    keep_score <- colSums(!is.finite(scores)) == 0
    keep_loading <- colSums(!is.finite(loadings)) == 0
    keep <- keep_score & keep_loading
    
    scores <- scores[, keep, drop = FALSE]
    loadings <- loadings[, keep, drop = FALSE]
    
    k <- min(ncol(scores), ncol(loadings))
    if (k <= 0) stop("No valid PCs remained after filtering.")
    
    list(
      scores   = scores[, seq_len(k), drop = FALSE],
      loadings = loadings[, seq_len(k), drop = FALSE],
      cutoff = if (is.null(res$cutoff)) NULL else res$cutoff
    )
  }, error = function(e) {
    message(sprintf("Method %s failed: %s", method_name, e$message))
    NULL
  })
}

# -----------------------------
# Clustering
# -----------------------------
run_kmeans_once <- function(emb, k, seed = 1L) {
  if (is.null(emb)) return(NULL)
  set.seed(seed)
  fit <- kmeans(emb, centers = k, nstart = 20)
  list(cluster = as.integer(fit$cluster), k = as.integer(k))
}

run_louvain_once <- function(emb, knn = 15L, seed = 1L) {
  if (is.null(emb)) return(NULL)
  if (nrow(emb) < 3) return(NULL)
  
  set.seed(seed)
  
  emb2 <- scale(emb)
  emb2[!is.finite(emb2)] <- 0
  
  k_use <- max(2L, min(knn, nrow(emb2) - 1L))
  # scran expects features x observations unless transposed=TRUE.
  g <- scran::buildSNNGraph(t(emb2), k = k_use, type = "jaccard")
  lv <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
  cl <- as.integer(lv$membership)
  
  list(cluster = cl, k = as.integer(length(unique(cl))))
}

run_gmm_once <- function(emb, k_max = 15L, seed = 1L) {
  if (is.null(emb)) return(NULL)
  if (nrow(emb) < 3) return(NULL)
  
  set.seed(seed)
  
  G_use <- 1:min(k_max, nrow(emb) - 1L)
  fit <- tryCatch(
    mclust::Mclust(emb, G = G_use, verbose = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit) || is.null(fit$classification)) return(NULL)
  
  list(
    cluster = as.integer(fit$classification),
    k = as.integer(fit$G)
  )
}

# Repeat clustering multiple times on the same embedding
get_clusterings_multi <- function(scores20, true_k, seed_base, cluster_reps = 10L) {
  emb10 <- select_pc(scores20, 10)
  emb20 <- select_pc(scores20, 20)
  
  reps_out <- vector("list", cluster_reps)
  time_out <- rep(NA_real_, cluster_reps)
  
  for (cr in seq_len(cluster_reps)) {
    t0 <- Sys.time()
    sb <- seed_base + cr * 1000L
    
    tmp <- list()
    
    if (RUN_KMEANS) {
      tmp$kmeans_pc10 <- tryCatch(
        run_kmeans_once(emb10, k = true_k, seed = sb + 10L),
        error = function(e) NULL
      )
      tmp$kmeans_pc20 <- tryCatch(
        run_kmeans_once(emb20, k = true_k, seed = sb + 20L),
        error = function(e) NULL
      )
    } else {
      tmp$kmeans_pc10 <- NULL
      tmp$kmeans_pc20 <- NULL
    }
    
    if (RUN_LOUVAIN) {
      tmp$louvain_pc10 <- tryCatch(
        run_louvain_once(emb10, knn = LOUVAIN_K, seed = sb + 30L),
        error = function(e) NULL
      )
      tmp$louvain_pc20 <- tryCatch(
        run_louvain_once(emb20, knn = LOUVAIN_K, seed = sb + 40L),
        error = function(e) NULL
      )
    } else {
      tmp$louvain_pc10 <- NULL
      tmp$louvain_pc20 <- NULL
    }
    
    if (RUN_GMM) {
      tmp$gmm_pc10 <- tryCatch(
        run_gmm_once(emb10, k_max = KMAX_MODEL, seed = sb + 50L),
        error = function(e) NULL
      )
      tmp$gmm_pc20 <- tryCatch(
        run_gmm_once(emb20, k_max = KMAX_MODEL, seed = sb + 60L),
        error = function(e) NULL
      )
    } else {
      tmp$gmm_pc10 <- NULL
      tmp$gmm_pc20 <- NULL
    }
    
    reps_out[[cr]] <- tmp
    time_out[cr] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }
  
  list(
    reps = reps_out,
    time = time_out
  )
}

# -----------------------------
# Metric block
# -----------------------------
algo_metric_block <- function(ref_obj, sub_obj, idx_sub, truth_full, truth_sub) {
  out <- list(
    ARI_truth_fullall   = NA_real_,
    NMI_truth_fullall   = NA_real_,
    ARI_truth_fullsub   = NA_real_,
    NMI_truth_fullsub   = NA_real_,
    ARI_truth_subset    = NA_real_,
    NMI_truth_subset    = NA_real_,
    ARI_full_vs_subset  = NA_real_,
    NMI_full_vs_subset  = NA_real_,
    k_full              = NA_integer_,
    k_subset            = NA_integer_
  )
  
  if (!is.null(ref_obj)) {
    out$ARI_truth_fullall <- ari_safe(ref_obj$cluster, truth_full)
    out$NMI_truth_fullall <- nmi_safe(ref_obj$cluster, truth_full)
    out$k_full <- ref_obj$k
  }
  
  if (!is.null(sub_obj)) {
    out$ARI_truth_subset <- ari_safe(sub_obj$cluster, truth_sub)
    out$NMI_truth_subset <- nmi_safe(sub_obj$cluster, truth_sub)
    out$k_subset <- sub_obj$k
  }
  
  if (!is.null(ref_obj)) {
    ref_sub <- ref_obj$cluster[idx_sub]
    out$ARI_truth_fullsub <- ari_safe(ref_sub, truth_sub)
    out$NMI_truth_fullsub <- nmi_safe(ref_sub, truth_sub)
    
    if (!is.null(sub_obj)) {
      out$ARI_full_vs_subset <- ari_safe(ref_sub, sub_obj$cluster)
      out$NMI_full_vs_subset <- nmi_safe(ref_sub, sub_obj$cluster)
    }
  }
  
  out
}

pc_similarity_block <- function(load_full, load_sub) {
  out <- list()
  
  pcwise20 <- pcwise_abs_cosine(load_full, load_sub, max_pc = 20L)
  out <- c(out, as.list(pcwise20))
  
  sim10 <- as.numeric(pcwise20[seq_len(10)])
  out$sim_mean_pc10 <- mean(sim10, na.rm = TRUE)
  out$sim_min_pc10  <- min(sim10, na.rm = TRUE)
  
  sim20 <- as.numeric(pcwise20[seq_len(20)])
  out$sim_mean_pc20 <- mean(sim20, na.rm = TRUE)
  out$sim_min_pc20  <- min(sim20, na.rm = TRUE)
  
  sub10 <- subspace_similarity(load_full, load_sub, p = 10L)
  sub20 <- subspace_similarity(load_full, load_sub, p = 20L)
  
  out <- c(out, append_subspace_prefix(sub10, "subspace_pc10_"))
  out <- c(out, append_subspace_prefix(sub20, "subspace_pc20_"))
  
  out
}

# -----------------------------
# Main
# -----------------------------
if (!file.exists(CSV_PATH)) stop("Pancreas data file not found: ", CSV_PATH)
if (dir.exists(OUT_DIR) && length(list.files(OUT_DIR, all.files = TRUE, no.. = TRUE))) {
  stop("Refusing to overwrite non-empty output directory: ", OUT_DIR)
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
write.csv(data.frame(
  code_version = "20260729_public_pcp_xv", dataset = "Pancreas", hvg = TOP_HVG, subset_n = SUBSET_N, repeats = REPS,
  cluster_repeats = CLUSTER_REPS, methods = paste(METHOD_NAMES, collapse = ","),
  calibration = "pairwise_empirical_quantile", alpha = PAIRWISE_ALPHA,
  base_seed = BASE_SEED, workers = WORKERS, data_root = DATA_ROOT,
  pcp_loading_source = "PCA_on_PCP_low_rank_L",
  pcp_score_construction = "X_times_V",
  stringsAsFactors = FALSE
), file.path(OUT_DIR, "run_manifest.csv"), row.names = FALSE)
PARTIAL_DIR <- file.path(OUT_DIR, "partial_reps")
dir.create(PARTIAL_DIR, recursive = TRUE, showWarnings = FALSE)

set.seed(BASE_SEED)
methods_all <- build_methods(MAX_PC)
methods <- methods_all[METHOD_NAMES]

cat("Reading real data...\n")
obj <- read_real_csv_pancreas(CSV_PATH, LABEL_COL)

counts_raw <- obj$counts
labels_raw <- obj$labels
labels_int <- to_int(labels_raw)

cat(sprintf("Raw data: %d cells x %d genes\n", nrow(counts_raw), ncol(counts_raw)))

cat("Preprocessing full data...\n")
prep <- preprocess_real_counts(counts_raw, top_hvg = TOP_HVG)
X_full <- prep$X

if (nrow(X_full) != length(labels_int)) {
  stop("Number of cells in X_full does not match label length.")
}
if (SUBSET_N > nrow(X_full)) {
  stop(sprintf("SUBSET_N=%d is larger than total cells=%d", SUBSET_N, nrow(X_full)))
}

k_true <- length(unique(labels_int))
cat(sprintf("After preprocessing: %d cells x %d genes | true k = %d\n",
            nrow(X_full), ncol(X_full), k_true))
cat(sprintf("[HVG] top variable genes selected: %d\n", ncol(X_full)))
cat(sprintf("[CLUSTER_REPS] repeated clustering per FULL/subset = %d\n", CLUSTER_REPS))

# -----------------------------
# Parallel plan
# -----------------------------
plan(multisession, workers = WORKERS)
cat(sprintf("Using workers: %d | nbrOfWorkers(): %d\n", WORKERS, nbrOfWorkers()))
print(plan())

# -----------------------------
# FULL reference: embedding + repeated clustering
# -----------------------------
cat("Computing FULL reference embeddings and repeated clusterings...\n")

full_ref <- list()

for (m in names(methods)) {
  cat(sprintf("  [FULL] %s ...\n", m))
  
  t0_embed <- Sys.time()
  fit_full <- safe_method_run(methods[[m]], X_full, method_name = m)
  dt_embed <- as.numeric(difftime(Sys.time(), t0_embed, units = "secs"))
  
  clu_full <- NULL
  if (!is.null(fit_full)) {
    clu_full <- get_clusterings_multi(
      scores20 = fit_full$scores,
      true_k = k_true,
      seed_base = BASE_SEED + 100000L + match(m, SUPPORTED_METHODS) * 1000L,
      cluster_reps = CLUSTER_REPS
    )
  } else {
    clu_full <- list(
      reps = replicate(CLUSTER_REPS, list(
        kmeans_pc10 = NULL, kmeans_pc20 = NULL,
        louvain_pc10 = NULL, louvain_pc20 = NULL,
        gmm_pc10 = NULL, gmm_pc20 = NULL
      ), simplify = FALSE),
      time = rep(NA_real_, CLUSTER_REPS)
    )
  }
  
  full_ref[[m]] <- list(
    scores = if (is.null(fit_full)) NULL else fit_full$scores,
    loadings = if (is.null(fit_full)) NULL else fit_full$loadings,
    cutoff = if (is.null(fit_full)) NULL else fit_full$cutoff,
    clustering_reps = clu_full$reps,
    time_full_embed = dt_embed,
    time_full_cluster = clu_full$time
  )
}

# -----------------------------
# Parallel sampling benchmark
# -----------------------------
cat("Running repeated subset benchmark in parallel...\n")

RNGkind("L'Ecuyer-CMRG")

rep_results <- future_lapply(
  seq_len(REPS),
  FUN = function(rep) {
    suppressPackageStartupMessages({
      library(mclust)
      library(aricode)
      library(rrcov)
      library(scran)
      library(igraph)
      library(rpca)
    })
    
    set.seed(BASE_SEED + rep)
    
    idx_sub <- sort(sample(seq_len(nrow(X_full)), SUBSET_N, replace = FALSE))
    X_sub <- X_full[idx_sub, , drop = FALSE]
    labels_sub <- labels_int[idx_sub]
    k_true_sub <- length(unique(labels_sub))
    
    rows_rep <- list()
    rr <- 0L
    
    for (m in names(methods)) {
      t0_embed <- Sys.time()
      fit_sub <- safe_method_run(methods[[m]], X_sub, method_name = paste0(m, "_subset"))
      dt_embed <- as.numeric(difftime(Sys.time(), t0_embed, units = "secs"))
      
      clu_sub <- NULL
      if (!is.null(fit_sub)) {
        clu_sub <- get_clusterings_multi(
          scores20 = fit_sub$scores,
          true_k = k_true_sub,
          seed_base = BASE_SEED + rep * 10000L + match(m, SUPPORTED_METHODS) * 100L,
          cluster_reps = CLUSTER_REPS
        )
      } else {
        clu_sub <- list(
          reps = replicate(CLUSTER_REPS, list(
            kmeans_pc10 = NULL, kmeans_pc20 = NULL,
            louvain_pc10 = NULL, louvain_pc20 = NULL,
            gmm_pc10 = NULL, gmm_pc20 = NULL
          ), simplify = FALSE),
          time = rep(NA_real_, CLUSTER_REPS)
        )
      }
      
      ref_m <- full_ref[[m]]
      
      for (cr in seq_len(CLUSTER_REPS)) {
        ref_clu <- ref_m$clustering_reps[[cr]]
        sub_clu <- clu_sub$reps[[cr]]
        
        row_out <- list(
          Dataset      = "SAMPLE",
          Replicate    = rep,          # outer subset replicate
          ClusterRep   = cr,           # inner clustering replicate
          SampleSize   = SUBSET_N,
          FullN        = nrow(X_full),
          GeneN        = ncol(X_full),
          TrueK        = k_true_sub,
          TrueK_full   = k_true,
          TrueK_subset = k_true_sub,
          Method       = m,
          Calibration  = if (m %in% c("Winsor", "Quad", "Ball", "Shell", "LR")) "pairwise_empirical_quantile" else "none",
          Alpha        = if (m %in% c("Winsor", "Quad", "Ball", "Shell", "LR")) PAIRWISE_ALPHA else NA_real_,
          Time_full_embed   = ref_m$time_full_embed,
          Time_full_cluster = ref_m$time_full_cluster[cr],
          Time_embed        = dt_embed,
          Time_cluster      = clu_sub$time[cr],
          Time_total        = dt_embed + clu_sub$time[cr]
        )
        row_out <- c(row_out, cutoff_columns(ref_m$cutoff, "full_"),
                     cutoff_columns(if (is.null(fit_sub)) NULL else fit_sub$cutoff, "subset_"))
        
        row_out <- c(row_out, append_prefix(
          algo_metric_block(
            ref_obj    = ref_clu$kmeans_pc10,
            sub_obj    = sub_clu$kmeans_pc10,
            idx_sub    = idx_sub,
            truth_full = labels_int,
            truth_sub  = labels_sub
          ),
          "kmeans_pc10_"
        ))
        
        row_out <- c(row_out, append_prefix(
          algo_metric_block(
            ref_obj    = ref_clu$kmeans_pc20,
            sub_obj    = sub_clu$kmeans_pc20,
            idx_sub    = idx_sub,
            truth_full = labels_int,
            truth_sub  = labels_sub
          ),
          "kmeans_pc20_"
        ))
        
        row_out <- c(row_out, append_prefix(
          algo_metric_block(
            ref_obj    = ref_clu$louvain_pc10,
            sub_obj    = sub_clu$louvain_pc10,
            idx_sub    = idx_sub,
            truth_full = labels_int,
            truth_sub  = labels_sub
          ),
          "louvain_pc10_"
        ))
        
        row_out <- c(row_out, append_prefix(
          algo_metric_block(
            ref_obj    = ref_clu$louvain_pc20,
            sub_obj    = sub_clu$louvain_pc20,
            idx_sub    = idx_sub,
            truth_full = labels_int,
            truth_sub  = labels_sub
          ),
          "louvain_pc20_"
        ))
        
        row_out <- c(row_out, append_prefix(
          algo_metric_block(
            ref_obj    = ref_clu$gmm_pc10,
            sub_obj    = sub_clu$gmm_pc10,
            idx_sub    = idx_sub,
            truth_full = labels_int,
            truth_sub  = labels_sub
          ),
          "gmm_pc10_"
        ))
        
        row_out <- c(row_out, append_prefix(
          algo_metric_block(
            ref_obj    = ref_clu$gmm_pc20,
            sub_obj    = sub_clu$gmm_pc20,
            idx_sub    = idx_sub,
            truth_full = labels_int,
            truth_sub  = labels_sub
          ),
          "gmm_pc20_"
        ))
        
        row_out <- c(
          row_out,
          pc_similarity_block(
            load_full = ref_m$loadings,
            load_sub  = if (is.null(fit_sub)) NULL else fit_sub$loadings
          )
        )
        
        rr <- rr + 1L
        rows_rep[[rr]] <- as.data.frame(row_out, stringsAsFactors = FALSE)
      }
    }
    
    rep_df <- do.call(rbind, rows_rep)
    write.csv(rep_df, file.path(PARTIAL_DIR, sprintf("rep_%02d.csv", rep)), row.names = FALSE)
    rep_df
  },
  future.seed = TRUE,
  future.packages = c("mclust", "aricode", "rrcov", "scran", "igraph", "rpca"),
  future.globals = c(
    "X_full", "labels_int", "full_ref", "methods", "PARTIAL_DIR",
    "SUBSET_N", "BASE_SEED", "REPS", "CLUSTER_REPS", "k_true", "PAIRWISE_ALPHA", "SUPPORTED_METHODS",
    "RUN_KMEANS", "RUN_LOUVAIN", "RUN_GMM",
    "LOUVAIN_K", "KMAX_MODEL",
    "safe_method_run", "get_clusterings_multi",
    "algo_metric_block", "append_prefix", "cutoff_columns",
    "pc_similarity_block", "ari_safe", "nmi_safe",
    "select_pc", "select_loading", "normalize_columns",
    "pcwise_abs_cosine", "subspace_similarity",
    "append_subspace_prefix",
    "run_kmeans_once", "run_louvain_once", "run_gmm_once",
    "build_methods", "pcp_decompose_fast", "pca_grid_fit_safe",
    "function_winsor_vectorized", "function_quad_vectorized",
    "function_ball_vectorized", "function_shell_vectorized",
    "function_lr_vectorized", "function_tau_vectorized", "pairwise_weighted_eigen"
  )
)

df_sample <- do.call(rbind, rep_results)

# -----------------------------
# FULL baseline table
# -----------------------------
full_rows <- list()
fi <- 0L

for (m in names(methods)) {
  ref_m <- full_ref[[m]]
  
  get_metric_full <- function(obj, truth) {
    if (is.null(obj)) {
      return(list(ARI = NA_real_, NMI = NA_real_, k = NA_integer_))
    }
    list(
      ARI = ari_safe(obj$cluster, truth),
      NMI = nmi_safe(obj$cluster, truth),
      k   = obj$k
    )
  }
  
  for (cr in seq_len(CLUSTER_REPS)) {
    fi <- fi + 1L
    clu <- ref_m$clustering_reps[[cr]]
    
    km10 <- get_metric_full(clu$kmeans_pc10, labels_int)
    km20 <- get_metric_full(clu$kmeans_pc20, labels_int)
    lv10 <- get_metric_full(clu$louvain_pc10, labels_int)
    lv20 <- get_metric_full(clu$louvain_pc20, labels_int)
    gm10 <- get_metric_full(clu$gmm_pc10, labels_int)
    gm20 <- get_metric_full(clu$gmm_pc20, labels_int)
    
    full_rows[[fi]] <- data.frame(
      Dataset = "FULL",
      Replicate = NA_integer_,
      ClusterRep = cr,
      SampleSize = nrow(X_full),
      FullN = nrow(X_full),
      GeneN = ncol(X_full),
      TrueK = k_true,
      TrueK_full = k_true,
      TrueK_subset = NA_integer_,
      Method = m,
      Calibration = if (m %in% c("Winsor", "Quad", "Ball", "Shell", "LR")) "pairwise_empirical_quantile" else "none",
      Alpha = if (m %in% c("Winsor", "Quad", "Ball", "Shell", "LR")) PAIRWISE_ALPHA else NA_real_,
      Time_full_embed = ref_m$time_full_embed,
      Time_full_cluster = ref_m$time_full_cluster[cr],
      Time = ref_m$time_full_embed + ref_m$time_full_cluster[cr],
      
      ARI_kmeans_pc10 = km10$ARI,
      NMI_kmeans_pc10 = km10$NMI,
      k_kmeans_pc10   = km10$k,
      
      ARI_kmeans_pc20 = km20$ARI,
      NMI_kmeans_pc20 = km20$NMI,
      k_kmeans_pc20   = km20$k,
      
      ARI_louvain_pc10 = lv10$ARI,
      NMI_louvain_pc10 = lv10$NMI,
      k_louvain_pc10   = lv10$k,
      
      ARI_louvain_pc20 = lv20$ARI,
      NMI_louvain_pc20 = lv20$NMI,
      k_louvain_pc20   = lv20$k,
      
      ARI_gmm_pc10 = gm10$ARI,
      NMI_gmm_pc10 = gm10$NMI,
      k_gmm_pc10   = gm10$k,
      
      ARI_gmm_pc20 = gm20$ARI,
      NMI_gmm_pc20 = gm20$NMI,
      k_gmm_pc20   = gm20$k,
      
      stringsAsFactors = FALSE
    )
  }
}

df_full <- do.call(rbind, full_rows)

# -----------------------------
# tidy sample-level ARI/NMI columns
# -----------------------------
rename_sample_metrics <- function(df) {
  df$ARI_kmeans_pc10 <- df$kmeans_pc10_ARI_truth_subset
  df$NMI_kmeans_pc10 <- df$kmeans_pc10_NMI_truth_subset
  df$k_kmeans_pc10   <- df$kmeans_pc10_k_subset
  
  df$ARI_kmeans_pc20 <- df$kmeans_pc20_ARI_truth_subset
  df$NMI_kmeans_pc20 <- df$kmeans_pc20_NMI_truth_subset
  df$k_kmeans_pc20   <- df$kmeans_pc20_k_subset
  
  df$ARI_louvain_pc10 <- df$louvain_pc10_ARI_truth_subset
  df$NMI_louvain_pc10 <- df$louvain_pc10_NMI_truth_subset
  df$k_louvain_pc10   <- df$louvain_pc10_k_subset
  
  df$ARI_louvain_pc20 <- df$louvain_pc20_ARI_truth_subset
  df$NMI_louvain_pc20 <- df$louvain_pc20_NMI_truth_subset
  df$k_louvain_pc20   <- df$louvain_pc20_k_subset
  
  df$ARI_gmm_pc10 <- df$gmm_pc10_ARI_truth_subset
  df$NMI_gmm_pc10 <- df$gmm_pc10_NMI_truth_subset
  df$k_gmm_pc10   <- df$gmm_pc10_k_subset
  
  df$ARI_gmm_pc20 <- df$gmm_pc20_ARI_truth_subset
  df$NMI_gmm_pc20 <- df$gmm_pc20_NMI_truth_subset
  df$k_gmm_pc20   <- df$gmm_pc20_k_subset
  
  df
}

df_sample <- rename_sample_metrics(df_sample)

# -----------------------------
# Optional summary tables
# -----------------------------
summary_mean_by_method <- function(df, value_cols) {
  sp <- split(df, df$Method)
  out <- lapply(names(sp), function(m) {
    sub <- sp[[m]]
    ans <- data.frame(Method = m, stringsAsFactors = FALSE)
    for (cc in value_cols) {
      ans[[paste0(cc, "_mean")]] <- mean(sub[[cc]], na.rm = TRUE)
      ans[[paste0(cc, "_median")]] <- median(sub[[cc]], na.rm = TRUE)
      ans[[paste0(cc, "_sd")]] <- sd(sub[[cc]], na.rm = TRUE)
    }
    ans
  })
  do.call(rbind, out)
}

# -----------------------------
# save outputs
# -----------------------------
FULL_CSV <- file.path(OUT_DIR, "real_full_results_repeated_clustering.csv")
SAMP_CSV <- file.path(OUT_DIR, "real_sampling_results_repeated_clustering.csv")

FULL_SUMMARY_CSV <- file.path(OUT_DIR, "real_full_results_repeated_clustering_summary.csv")
SAMP_SUMMARY_CSV <- file.path(OUT_DIR, "real_sampling_results_repeated_clustering_summary.csv")

write.csv(df_full, FULL_CSV, row.names = FALSE)
write.csv(df_sample, SAMP_CSV, row.names = FALSE)
write.csv(df_sample, OUT_CSV, row.names = FALSE)

full_summary <- summary_mean_by_method(
  df_full,
  value_cols = c("ARI_gmm_pc10", "ARI_gmm_pc20",
                 "ARI_kmeans_pc10", "ARI_kmeans_pc20",
                 "ARI_louvain_pc10", "ARI_louvain_pc20",
                 "NMI_gmm_pc10", "NMI_gmm_pc20",
                 "NMI_kmeans_pc10", "NMI_kmeans_pc20",
                 "NMI_louvain_pc10", "NMI_louvain_pc20")
)

sample_summary <- summary_mean_by_method(
  df_sample,
  value_cols = c("ARI_gmm_pc10", "ARI_gmm_pc20",
                 "ARI_kmeans_pc10", "ARI_kmeans_pc20",
                 "ARI_louvain_pc10", "ARI_louvain_pc20",
                 "NMI_gmm_pc10", "NMI_gmm_pc20",
                 "NMI_kmeans_pc10", "NMI_kmeans_pc20",
                 "NMI_louvain_pc10", "NMI_louvain_pc20",
                 "gmm_pc10_ARI_full_vs_subset", "gmm_pc20_ARI_full_vs_subset",
                 "kmeans_pc10_ARI_full_vs_subset", "kmeans_pc20_ARI_full_vs_subset",
                 "louvain_pc10_ARI_full_vs_subset", "louvain_pc20_ARI_full_vs_subset")
)

write.csv(full_summary, FULL_SUMMARY_CSV, row.names = FALSE)
write.csv(sample_summary, SAMP_SUMMARY_CSV, row.names = FALSE)

cat(sprintf(
  "✅ Done.\nFULL results: %s\nSAMPLE results: %s\nFULL summary: %s\nSAMPLE summary: %s\n",
  FULL_CSV, SAMP_CSV, FULL_SUMMARY_CSV, SAMP_SUMMARY_CSV
))
