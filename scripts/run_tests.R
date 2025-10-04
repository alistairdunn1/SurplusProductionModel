# Portable test runner for SurplusProductionModel
suppressPackageStartupMessages({
  if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("devtools is required to run tests. Please install it: install.packages('devtools')")
  }
})

# Determine the directory of this script, then set package_dir to its parent
get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grepl("^--file=", args)]
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE))
  }
  # Fallback to current working directory when running interactively
  return(normalizePath(getwd(), winslash = "/", mustWork = TRUE))
}

script_path <- get_script_path()
script_dir <- dirname(script_path)
package_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
message("Package dir: ", package_dir)

setwd(package_dir)

devtools::load_all(quiet = TRUE)
message("Loaded package. Running tests...")

if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("testthat is required to run tests. Please install it: install.packages('testthat')")
}

# Prefer devtools::test() which uses testthat under the hood
res <- tryCatch({
  devtools::test(stop_on_failure = FALSE, stop_on_warning = FALSE)
}, error = function(e) e)

if (inherits(res, "error")) {
  message("Tests failed with error: ", conditionMessage(res))
  q(status = 1)
} else {
  message("Tests completed. Review output above.")
}
