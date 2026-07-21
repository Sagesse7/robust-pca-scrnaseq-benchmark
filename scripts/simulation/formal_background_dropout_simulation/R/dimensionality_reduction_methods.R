# Shared dimensionality-reduction methods used by Simulation 1 and Simulation 2.
#
# Proposed generalized Kendall methods use pairwise empirical-quantile
# calibration. The squared weights returned here are the quantities that enter
# K_hat = 1/[n(n-1)] sum_{i != j} w_ij (Xi-Xj)(Xi-Xj)'.

pairwise_geometry <- function(X) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  sq <- rowSums(X^2)
  D2 <- pmax(outer(sq, sq, "+") - 2 * tcrossprod(X), 0)
  diag(D2) <- 0
  D <- sqrt(D2)
  list(D = D, D2 = D2, upper = D[upper.tri(D)])
}

pairwise_empirical_cutoffs <- function(geometry, alpha = 0.05) {
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be strictly between 0 and 1.")
  }
  q <- function(p) as.numeric(stats::quantile(
    geometry$upper, probs = p, names = FALSE, type = 7, na.rm = TRUE
  ))
  list(
    calibration = "pairwise_empirical_quantile",
    alpha = alpha,
    Q1 = q(alpha / 2),
    Q2 = q(1 - alpha),
    Q3 = q(1 - alpha / 2),
    Q3star = q(1 - alpha / 4)
  )
}

squared_pairwise_weights <- function(method, geometry, cutoff = NULL) {
  D <- geometry$D
  D2 <- geometry$D2

  if (method == "Tau") {
    W <- matrix(0, nrow(D), ncol(D))
    W[D2 > 0] <- 1 / D2[D2 > 0]
  } else if (method == "Winsor") {
    W <- ifelse(D <= cutoff$Q2, 1, cutoff$Q2^2 / D2)
  } else if (method == "Quad") {
    W <- ifelse(D <= cutoff$Q2, 1, cutoff$Q2^4 / D2^2)
  } else if (method == "Ball") {
    W <- ifelse(D <= cutoff$Q2, 1, 0)
  } else if (method == "Shell") {
    W <- ifelse(D >= cutoff$Q1 & D <= cutoff$Q3, 1, 0)
  } else if (method == "LR") {
    denom <- cutoff$Q3star - cutoff$Q2
    if (!is.finite(denom) || denom <= 0) stop("LR requires Q3star > Q2.")
    W <- ifelse(
      D <= cutoff$Q2,
      1,
      ifelse(D <= cutoff$Q3star, ((cutoff$Q3star - D) / denom)^2, 0)
    )
  } else {
    stop("Unsupported generalized Kendall method: ", method)
  }

  W[!is.finite(W)] <- 0
  diag(W) <- 0
  (W + t(W)) / 2
}

weighted_pairwise_matrix <- function(X, W) {
  n <- nrow(X)
  S <- rowSums(W)
  K <- 2 * (crossprod(X, S * X) - crossprod(X, W %*% X)) / (n * (n - 1))
  (K + t(K)) / 2
}

top_symmetric_eigen <- function(K, k = 20L) {
  k <- min(as.integer(k), ncol(K))
  if (k < 1L) stop("k must be positive.")

  if (requireNamespace("RSpectra", quietly = TRUE) && k < ncol(K) - 1L) {
    fit <- RSpectra::eigs_sym(K, k = k, which = "LA")
    ord <- order(fit$values, decreasing = TRUE)
    values <- fit$values[ord]
    vectors <- fit$vectors[, ord, drop = FALSE]
    backend <- "RSpectra::eigs_sym"
  } else {
    fit <- eigen(K, symmetric = TRUE)
    values <- fit$values[seq_len(k)]
    vectors <- fit$vectors[, seq_len(k), drop = FALSE]
    backend <- "base::eigen"
  }

  if (!is.numeric(vectors) || any(!is.finite(vectors))) {
    stop("Invalid eigenvectors returned.")
  }
  list(values = values, loadings = vectors, backend = backend)
}

fit_generalized_kendall <- function(
    X,
    method,
    k = 20L,
    alpha = 0.05,
    geometry = NULL,
    return_weights = FALSE
) {
  method <- match.arg(method, c("Tau", "Winsor", "Quad", "Ball", "Shell", "LR"))
  if (is.null(geometry)) geometry <- pairwise_geometry(X)
  cutoff <- if (method == "Tau") NULL else pairwise_empirical_cutoffs(geometry, alpha)
  W <- squared_pairwise_weights(method, geometry, cutoff)
  fit <- top_symmetric_eigen(weighted_pairwise_matrix(X, W), k)
  fit$method <- method
  fit$calibration <- if (method == "Tau") "none" else cutoff$calibration
  fit$alpha <- if (method == "Tau") NA_real_ else alpha
  fit$cutoff <- cutoff
  if (return_weights) fit$weights <- W
  fit
}

fit_pca <- function(X, k = 20L) {
  fit <- prcomp(X, center = FALSE, scale. = FALSE)
  k <- min(as.integer(k), ncol(fit$rotation))
  list(
    method = "PCA", calibration = "none", alpha = NA_real_,
    loadings = fit$rotation[, seq_len(k), drop = FALSE],
    values = fit$sdev[seq_len(k)]^2, backend = "stats::prcomp"
  )
}

fit_hubert <- function(X, k = 20L) {
  if (!requireNamespace("rrcov", quietly = TRUE)) stop("Package rrcov is required.")
  k <- min(as.integer(k), ncol(X), nrow(X) - 1L)
  fit <- rrcov::PcaHubert(X, k = k, kmax = k, scale = FALSE)
  list(
    method = "Hubert", calibration = "none", alpha = NA_real_,
    loadings = as.matrix(fit@loadings)[, seq_len(k), drop = FALSE],
    values = as.numeric(fit@eigenvalues)[seq_len(k)], backend = "rrcov::PcaHubert"
  )
}

fit_grid <- function(X, k = 20L) {
  if (!requireNamespace("rrcov", quietly = TRUE)) stop("Package rrcov is required.")
  p <- ncol(X)
  k <- min(as.integer(k), p, nrow(X) - 1L)
  fit <- rrcov::PcaGrid(
    as.data.frame(X), k = k, kmax = k,
    center = rep(0, p), scale = rep(1, p), method = "mad"
  )
  loadings <- as.matrix(fit@loadings)[, seq_len(k), drop = FALSE]
  list(
    method = "Grid", calibration = "none", alpha = NA_real_,
    loadings = loadings,
    values = as.numeric(fit@eigenvalues)[seq_len(k)], backend = "rrcov::PcaGrid",
    grid_scale_method = "mad",
    grid_computed_k = as.integer(fit@k),
    grid_orthogonality_error = max(abs(crossprod(loadings) - diag(ncol(loadings))))
  )
}

fit_pcp <- function(X, k = 20L, lambda = 1 / sqrt(max(dim(X)))) {
  if (!requireNamespace("rpca", quietly = TRUE)) stop("Package rpca is required.")
  rp <- rpca::rpca(
    X, lambda = lambda, term.delta = 1e-5, max.iter = 500, trace = FALSE
  )
  fit <- prcomp(rp$L, center = FALSE, scale. = FALSE)
  k <- min(as.integer(k), ncol(fit$rotation))
  relative_reconstruction_error <- sqrt(sum((X - rp$L - rp$S)^2)) / sqrt(sum(X^2))
  list(
    method = "PCP", calibration = "none", alpha = NA_real_,
    loadings = fit$rotation[, seq_len(k), drop = FALSE],
    values = fit$sdev[seq_len(k)]^2, backend = "rpca::rpca+stats::prcomp",
    pcp_lambda = lambda,
    pcp_converged = isTRUE(rp$convergence$converged),
    pcp_iterations = as.integer(rp$convergence$iterations),
    pcp_final_delta = as.numeric(rp$convergence$final.delta),
    pcp_relative_reconstruction_error = relative_reconstruction_error,
    pcp_low_rank_estimate = length(rp$L.svd$d),
    pcp_sparse_fraction = mean(rp$S != 0)
  )
}

fit_diagnostic_fields <- function(fit) {
  wanted <- c(
    "grid_scale_method", "grid_computed_k", "grid_orthogonality_error",
    "pcp_lambda", "pcp_converged", "pcp_iterations", "pcp_final_delta",
    "pcp_relative_reconstruction_error", "pcp_low_rank_estimate",
    "pcp_sparse_fraction"
  )
  fit[intersect(wanted, names(fit))]
}

all_method_names <- function() {
  c("Grid", "Hubert", "PCP", "PCA", "Tau", "Winsor", "Quad", "Ball", "Shell", "LR")
}

proposed_method_names <- function() c("Winsor", "Quad", "Ball", "Shell", "LR")

fit_dimensionality_reduction <- function(
    X,
    method,
    k = 20L,
    alpha = 0.05,
    geometry = NULL,
    return_weights = FALSE
) {
  method <- match.arg(method, all_method_names())
  if (method == "PCA") return(fit_pca(X, k))
  if (method == "Hubert") return(fit_hubert(X, k))
  if (method == "Grid") return(fit_grid(X, k))
  if (method == "PCP") return(fit_pcp(X, k))
  fit_generalized_kendall(X, method, k, alpha, geometry, return_weights)
}

method_run_grid <- function(methods = all_method_names(), alphas = c(0.05, 0.10, 0.20)) {
  methods <- intersect(methods, all_method_names())
  if (!length(methods)) stop("No supported methods selected.")
  out <- lapply(methods, function(method) {
    if (method %in% proposed_method_names()) {
      data.frame(Method = method, Alpha = alphas)
    } else {
      data.frame(Method = method, Alpha = NA_real_)
    }
  })
  do.call(rbind, out)
}

safe_fit_dimensionality_reduction <- function(...) {
  tryCatch(
    {
      fit <- fit_dimensionality_reduction(...)
      if (!is.matrix(fit$loadings) || ncol(fit$loadings) < 1L ||
          any(!is.finite(fit$loadings)) || any(!is.finite(fit$values))) {
        stop("Method returned invalid loadings or eigenvalues.")
      }
      list(success = TRUE, fit = fit, failure_message = "")
    },
    error = function(e) {
      list(success = FALSE, fit = NULL, failure_message = conditionMessage(e))
    }
  )
}
