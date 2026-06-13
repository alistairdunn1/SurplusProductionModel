# Example script demonstrating Pella-Tomlinson model fitting

# Load the package in development mode
devtools::load_all(".")

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
cpue_data <- pmax(0.001, true_q * true_biomass * exp(rnorm(n_years, 0, 0.05))) # Reduced noise and minimum value

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

# Calculate reference points
cat("\nCalculating reference points...\n")
ref_points <- calculate_reference_points(fitted_model)
print(ref_points)

cat("\nExample completed successfully!\n")
