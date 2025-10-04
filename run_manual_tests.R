# Run testthat tests manually without full package loading
setwd("c:/Users/alist/OneDrive/Projects/Software/ATO_rTMB/SurplusProductionModel")

# Load required libraries
library(methods)
library(checkmate)
library(testthat)

cat("Loading S4 classes for testing...\n")

# Source files in order
source("R/generics.R")
source("R/classes.R")
source("R/constructors.R")

cat("Running testthat unit tests...\n\n")

# Run the specific test file
test_results <- test_file("tests/testthat/test-classes.R", reporter = "summary")

cat("\n")
cat("========================================\n")
cat("TEST SUMMARY\n")
cat("========================================\n")

if (length(test_results) > 0) {
  passed <- sum(sapply(test_results, function(x) x$nb))
  failed <- sum(sapply(test_results, function(x) length(x$failed)))

  cat("Tests passed:", passed, "\n")
  cat("Tests failed:", failed, "\n")

  if (failed == 0) {
    cat("🎉 All tests PASSED!\n")
  } else {
    cat("❌ Some tests FAILED!\n")
  }
} else {
  cat("⚠️  No test results returned\n")
}
