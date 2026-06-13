# Example script demonstrating Pella-Tomlinson model fitting

# Load the package
library(SurplusProductionModel)

# Create synthetic example data for demonstration
set.seed(123)
years <- 2010:2020
n_years <- length(years)

# "True" model parameters (for generating synthetic data)
true_r <- 0.3
true_K <- 5000
true_m <- 2.0 # Schaefer model
true_q <- 0.001
true_B_initial <- 4000

# Generate deterministic biomass trajectory
true_biomass <- numeric(n_years)
true_biomass[1] <- true_B_initial
catch_data <- rep(800, n_years) # Constant catch

for (t in 1:(n_years - 1)) {
  # Pella-Tomlinson production
  production <- true_r * true_biomass[t] * (1 - true_biomass[t] / true_K)
  true_biomass[t + 1] <- true_biomass[t] + production - catch_data[t]
}

# Generate CPUE with observation error
cpue_data <- true_q * true_biomass * exp(rnorm(n_years, 0, 0.1))

# Prepare data for model fitting
model_data <- list(
  cpue_data = data.frame(
    year = years,
    cpue = cpue_data
  ),
  catch_data = data.frame(
    year = years,
    catch = catch_data
  )
)

# Display the data
cat("Example CPUE and Catch Data:\n")
print(model_data$cpue_data)
print(model_data$catch_data)

# Fit the Pella-Tomlinson model
cat("\nFitting Pella-Tomlinson surplus production model...\n")

# Set options for faster convergence in example
options_list <- list(
  validate_data = TRUE,
  silent = TRUE,
  control = list(eval.max = 200, iter.max = 100)
)

# Fit the model
fitted_model <- fit_pella_tomlinson_model(
  data = model_data,
  options = options_list
)

# Display model results
cat("\nModel Fitting Results:\n")
print(fitted_model)

# Extract fitted parameters
fitted_params <- fitted_model$parameters
cat("\nFitted Parameters:\n")
cat("r (intrinsic growth rate):", round(fitted_params["r"], 4), "\n")
cat("K (carrying capacity):", round(fitted_params["K"], 1), "tonnes\n")
cat("m (shape parameter):", round(fitted_params["m"], 2), "\n")
cat("q (catchability):", round(fitted_params["q"], 6), "\n")

# Calculate reference points
cat("\nCalculating reference points...\n")
ref_points <- calculate_reference_points(fitted_model)
print(ref_points)

# Extract biomass estimates
cat("\nBiomass estimates:\n")
biomass_estimates <- estimate_biomass(fitted_model)
print(biomass_estimates)

# Compare with "true" values
cat("\nComparison with true values:\n")
cat("True r:", true_r, "vs Fitted r:", round(fitted_params["r"], 4), "\n")
cat("True K:", true_K, "vs Fitted K:", round(fitted_params["K"], 1), "\n")
cat("True m:", true_m, "vs Fitted m:", round(fitted_params["m"], 2), "\n")
cat("True q:", true_q, "vs Fitted q:", round(fitted_params["q"], 6), "\n")

# Model diagnostics
if ("likelihood" %in% names(fitted_model$results)) {
  cat("\nModel fit statistics:\n")
  cat("Negative log-likelihood:", round(fitted_model$results$likelihood, 2), "\n")
  cat("AIC:", round(fitted_model$results$aic, 2), "\n")
  cat("Convergence code:", fitted_model$results$convergence, "\n")
}

cat("\nExample completed successfully!\n")
