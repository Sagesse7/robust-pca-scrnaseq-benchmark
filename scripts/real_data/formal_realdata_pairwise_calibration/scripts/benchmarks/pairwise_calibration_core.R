# Pairwise empirical-quantile calibration for generalized Kendall PCA.

pairwise_weighted_eigen <- function(Y, method, alpha = 0.05) {
  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"
  n <- nrow(Y)
  if (n < 3L) stop("At least three observations are required.")
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must lie strictly between 0 and 1.")
  }

  sq_norms <- rowSums(Y^2)
  D2 <- outer(sq_norms, sq_norms, "+") - 2 * tcrossprod(Y)
  D2[D2 < 0 & D2 > -1e-8] <- 0
  if (any(D2 < 0, na.rm = TRUE)) stop("Negative squared distances exceeded tolerance.")
  D <- sqrt(D2)
  pair_index <- upper.tri(D)
  pairwise <- D[pair_index]
  pairwise <- pairwise[is.finite(pairwise)]
  if (!length(pairwise)) stop("No finite pairwise distances were available.")

  probs <- c(alpha / 2, 1 - alpha, 1 - alpha / 2, 1 - alpha / 4)
  q <- as.numeric(stats::quantile(pairwise, probs = probs, names = FALSE, type = 7))
  names(q) <- c("Q1", "Q2", "Q3", "Q3star")

  if (method == "Winsor") {
    W <- ifelse(D <= q[["Q2"]], 1, (q[["Q2"]] / D)^2)
  } else if (method == "Quad") {
    W <- ifelse(D <= q[["Q2"]], 1, (q[["Q2"]] / D)^4)
  } else if (method == "Ball") {
    W <- ifelse(D <= q[["Q2"]], 1, 0)
  } else if (method == "Shell") {
    W <- ifelse(D >= q[["Q1"]] & D <= q[["Q3"]], 1, 0)
  } else if (method == "LR") {
    width <- q[["Q3star"]] - q[["Q2"]]
    if (!is.finite(width) || width <= .Machine$double.eps) {
      W <- ifelse(D <= q[["Q2"]], 1, 0)
    } else {
      W <- ifelse(
        D <= q[["Q2"]], 1,
        ifelse(D <= q[["Q3star"]], ((q[["Q3star"]] - D) / width)^2, 0)
      )
    }
  } else {
    stop("Unsupported generalized method: ", method)
  }

  W[!is.finite(W)] <- 0
  diag(W) <- 0
  pair_weights <- W[pair_index]
  sum_w <- sum(pair_weights)
  sum_w2 <- sum(pair_weights^2)
  n_eff <- if (sum_w2 > 0) sum_w^2 / sum_w2 else 0

  row_weight <- rowSums(W)
  K <- 2 * (crossprod(Y, row_weight * Y) - crossprod(Y, W %*% Y)) /
    (n * (n - 1))
  K <- (K + t(K)) / 2
  if (any(!is.finite(K))) stop("The generalized Kendall matrix contains non-finite values.")

  eg <- eigen(K, symmetric = TRUE)
  eg$cutoff <- list(
    calibration = "pairwise_empirical_quantile",
    alpha = alpha,
    Q1 = unname(q[["Q1"]]),
    Q2 = unname(q[["Q2"]]),
    Q3 = unname(q[["Q3"]]),
    Q3star = unname(q[["Q3star"]]),
    mean_weight = mean(pair_weights),
    median_weight = stats::median(pair_weights),
    prop_weight_lt_1 = mean(pair_weights < 1 - 1e-12),
    prop_weight_lt_0.5 = mean(pair_weights < 0.5),
    prop_weight_zero = mean(pair_weights <= 1e-12),
    effective_pairs = n_eff,
    total_pairs = length(pair_weights)
  )
  eg
}

function_winsor_vectorized <- function(Y) {
  pairwise_weighted_eigen(Y, "Winsor", PAIRWISE_ALPHA)
}

function_quad_vectorized <- function(Y) {
  pairwise_weighted_eigen(Y, "Quad", PAIRWISE_ALPHA)
}

function_ball_vectorized <- function(Y) {
  pairwise_weighted_eigen(Y, "Ball", PAIRWISE_ALPHA)
}

function_shell_vectorized <- function(Y) {
  pairwise_weighted_eigen(Y, "Shell", PAIRWISE_ALPHA)
}

function_lr_vectorized <- function(Y) {
  pairwise_weighted_eigen(Y, "LR", PAIRWISE_ALPHA)
}

function_tau_vectorized <- function(Y) {
  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"
  n <- nrow(Y)
  if (n < 2L) stop("At least two observations are required.")

  sq_norms <- rowSums(Y^2)
  D2 <- outer(sq_norms, sq_norms, "+") - 2 * tcrossprod(Y)
  D2[D2 < 0 & D2 > -1e-8] <- 0
  if (any(D2 < 0, na.rm = TRUE)) stop("Negative squared distances exceeded tolerance.")

  W <- 1 / D2
  W[!is.finite(W)] <- 0
  diag(W) <- 0
  row_weight <- rowSums(W)
  K <- 2 * (crossprod(Y, row_weight * Y) - crossprod(Y, W %*% Y)) /
    (n * (n - 1))
  K <- (K + t(K)) / 2
  if (any(!is.finite(K))) stop("The Kendall matrix contains non-finite values.")
  eigen(K, symmetric = TRUE)
}
