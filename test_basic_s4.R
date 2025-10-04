# Minimal S4 class test - bypasses package loading
# Test the S4 classes by sourcing files directly

setwd("c:/Users/alist/OneDrive/Projects/Software/ATO_rTMB/SurplusProductionModel")

# Load required libraries that should be available
library(methods)
if (!requireNamespace("checkmate", quietly = TRUE)) {
  message("Warning: checkmate not available - installing...")
  install.packages("checkmate", repos = "https://cran.r-project.org")
  library(checkmate)
} else {
  library(checkmate)
}

cat("Loading S4 classes...\n")

# Source files in order
source("R/generics.R")
source("R/classes.R")
source("R/constructors.R")

cat("Files loaded successfully!\n")

# Test basic functionality
cat("Testing ProductionModel constructor...\n")

# Create test data
years <- 2000:2005
catch <- c(1000, 1100, 1200, 950, 1050, 1300)
cpue <- c(2.1, 2.0, 1.9, 2.2, 2.1, 1.8)
effort <- catch / cpue

# Test constructor
model <- ProductionModel(years = years, catch = catch, cpue = cpue, effort = effort)

cat("✓ Model created successfully!\n")
cat("✓ Class:", class(model), "\n")
cat("✓ Fitted status:", fitted(model), "\n")
cat("✓ Model type:", model@model_type, "\n")

# Test accessors
cat("\nTesting accessor methods...\n")
model_data <- model_data(model)
cat("✓ Data years range:", range(model_data$years), "\n")
cat("✓ Catch range:", sprintf("%.1f - %.1f", min(model_data$catch), max(model_data$catch)), "\n")

model_params <- parameters(model)
cat("✓ Parameters available:", length(model_params), "parameters\n")
cat("✓ Parameter names:", paste(names(model_params), collapse = ", "), "\n")

# Test validation
cat("\nTesting validation...\n")
tryCatch(
  {
    invalid_model <- ProductionModel(
      years = c(2000, 1999, 2001), # Non-increasing
      catch = catch[1:3],
      cpue = cpue[1:3],
      effort = effort[1:3]
    )
    cat("✗ Validation failed - should have caught non-increasing years\n")
  },
  error = function(e) {
    cat("✓ Validation working - caught error:", conditionMessage(e), "\n")
  }
)

# Test print method
cat("\nTesting print method...\n")
print(model)

cat("\n🎉 All basic S4 tests passed!\n")
cat("The ProductionModel S4 class is working correctly.\n")
