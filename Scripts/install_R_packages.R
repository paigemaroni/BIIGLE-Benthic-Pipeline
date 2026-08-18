#!/usr/bin/env Rscript

packages <- c(
  "readr", "dplyr", "tidyr", "ggplot2", "stringr", "vegan", "mgcv"
)

missing <- packages[
  !vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing) == 0) {
  cat("All required R packages are already installed.\n")
} else {
  cat("Installing missing packages:\n")
  cat(paste0("  - ", missing, collapse = "\n"), "\n\n")
  install.packages(missing)
}
