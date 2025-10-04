# Test script for S4 classes
setwd("c:/Users/alist/OneDrive/Projects/Software/ATO_rTMB/SurplusProductionModel")

# Load the package
devtools::load_all(".")

# Test basic functionality
cat("Testing ProductionModel constructor...\n")

# Create test data
years <- 2000:2010
catch <- c(1000, 1100, 1200, 950, 1050, 1300, 1150, 1250, 1400, 1350, 1200)
cpue <- c(2.1, 2.0, 1.9, 2.2, 2.1, 1.8, 2.0, 1.9, 1.7, 1.8, 1.9)
effort <- catch / cpue

# Test constructor
model <- ProductionModel(years = years, catch = catch, cpue = cpue, effort = effort)

cat("Model created successfully!\n")
cat("Class:", class(model), "\n")
cat("Fitted status:", fitted(model), "\n")
cat("Model type:", model@model_type, "\n")

# Test accessors
cat("\nTesting accessor methods...\n")
model_data <- data(model)
cat("Data years range:", range(model_data$years), "\n")
cat("Catch range:", range(model_data$catch), "\n")

model_params <- parameters(model)
cat("Parameters:", names(model_params), "\n")

# Test print method
cat("\nTesting print method...\n")
print(model)

cat("\nAll tests completed successfully!\n")
