args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[[1]] else "sessionInfo.txt"

packages <- c(
  "splatter", "SummarizedExperiment", "data.table", "mclust", "aricode",
  "scran", "igraph", "rrcov", "rpca", "RSpectra", "readxl", "digest"
)

lines <- c(
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("Platform: %s", R.version$platform),
  "",
  "Package versions:"
)
for (package in packages) {
  version <- if (requireNamespace(package, quietly = TRUE)) {
    as.character(utils::packageVersion(package))
  } else {
    "NOT INSTALLED"
  }
  lines <- c(lines, sprintf("%-24s %s", package, version))
}

lines <- c(lines, "", "Full sessionInfo():", capture.output(sessionInfo()))
writeLines(lines, output)
cat("Wrote reproducibility information to:", normalizePath(output), "\n")
