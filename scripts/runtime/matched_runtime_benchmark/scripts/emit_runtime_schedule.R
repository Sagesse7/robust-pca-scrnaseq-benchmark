#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) != 1L) stop("Supply exactly one --", name, "= value.")
  sub(paste0("^--", name, "="), "", hit)
}

replicate <- as.integer(arg_value("replicate"))
if (!is.finite(replicate) || replicate < 1L || replicate > 5L) {
  stop("replicate must be an integer from 1 to 5.")
}

methods <- c("PCA", "Grid", "Hubert", "PCP", "Tau",
             "Winsor", "Quad", "Ball", "Shell", "LR")
schedule <- expand.grid(
  SettingID = 1:7,
  Method = methods,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

# A deterministic, replicate-specific permutation prevents a fixed method or
# dimension from always being timed early or late in a long single-node job.
set.seed(270000L + replicate)
schedule <- schedule[sample.int(nrow(schedule)), , drop = FALSE]
write.table(schedule, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
