script_dir <- dirname(normalizePath(sub("^--file=", "", grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)[1])))
release_dir <- dirname(script_dir)
source(file.path(release_dir, "R", "dimensionality_reduction_methods.R"))

set.seed(260614)
X <- scale(matrix(rnorm(40 * 25), nrow = 40, ncol = 25))
geometry <- pairwise_geometry(X)
cutoff <- pairwise_empirical_cutoffs(geometry, alpha = 0.05)

stopifnot(
  cutoff$Q1 < cutoff$Q2,
  cutoff$Q2 < cutoff$Q3,
  cutoff$Q3 < cutoff$Q3star
)

for (method in c("Tau", proposed_method_names(), "PCA")) {
  result <- safe_fit_dimensionality_reduction(
    X, method = method, k = 5, alpha = 0.05, geometry = geometry,
    return_weights = method != "PCA"
  )
  stopifnot(result$success)
  stopifnot(identical(dim(result$fit$loadings), c(25L, 5L)))
  if (!is.null(result$fit$weights)) {
    stopifnot(max(abs(result$fit$weights - t(result$fit$weights))) < 1e-12)
    stopifnot(all(diag(result$fit$weights) == 0))
  }
}

if (requireNamespace("rrcov", quietly = TRUE)) {
  grid <- safe_fit_dimensionality_reduction(X, method = "Grid", k = 5)
  stopifnot(grid$success)
  stopifnot(grid$fit$grid_scale_method == "mad")
  stopifnot(grid$fit$grid_computed_k == 5L)
  stopifnot(grid$fit$grid_orthogonality_error < 1e-8)
}

if (requireNamespace("rpca", quietly = TRUE)) {
  pcp <- safe_fit_dimensionality_reduction(X, method = "PCP", k = 5)
  stopifnot(pcp$success)
  stopifnot(is.logical(pcp$fit$pcp_converged))
  stopifnot(is.finite(pcp$fit$pcp_relative_reconstruction_error))
  stopifnot(pcp$fit$pcp_sparse_fraction >= 0, pcp$fit$pcp_sparse_fraction <= 1)
}

cat("Method smoke test passed.\n")
