# SurplusProductionModel

## Overview

`SurplusProductionModel` is an R package that implements a spatial Pella-Tomlinson surplus production model for Antarctic toothfish (*Dissostichus mawsoni*) stock assessment. It supports an optional state-space formulation with process error for model fitting. This package serves as the operating model foundation for the the Surplus Production Model Management Strategy Evaluation (MSE) framework.

## Key Features

### Model Fitting

- **Pella-Tomlinson production function** with flexible shape parameter (*m*) (Pella and Tomlinson 1969)
  - Schaefer model (*m* = 2): symmetric production curve (Schaefer 1954)
  - Fox model (*m* = 1): asymmetric production curve with peak at lower biomass (Fox 1970)
- **Optional state-space framework** with process error in biomass dynamics when `options$process_noise = TRUE`, plus observation error in CPUE
- **RTMB integration** for automatic differentiation and fast optimization via `nlminb`
- **Spatial structure** supporting multiple management areas with optional movement; the movement kernel follows a gravity formulation (Turchin 1998)
- **Multi-index support** for multiple CPUE series per area

### Convergence Diagnostics

- **Multi-start optimization** (`n_starts` parameter) to find global minimum
- **Jitter test** (`jitter_test()`) to assess optimization reliability from perturbed starts
- **Retrospective analysis** (`retrospective_analysis()`) with Mohn's rho for systematic bias detection

### Uncertainty Quantification

- **Delta-method standard errors** via `sdreport` (Hessian-based)
- **Profile likelihood confidence intervals** (`profile_likelihood()`) for parameters and derived quantities
- **Bayesian MCMC sampling** (`bayesian_fit()`) via Stan's NUTS sampler using `tmbstan`
- **Posterior predictive checks** (`posterior_predictive_check()`) for model validation

### Diagnostics & Visualization

- **Reference points**: MSY, B`<sub>`MSY`</sub>`, F`<sub>`MSY`</sub>` with current B/B`<sub>`MSY`</sub>` ratio
- **Residual diagnostics** (`plot_residuals()`): QQ plots, residuals vs. fitted, histograms
- **Biomass plots** (`plot_biomass()`): trajectories with confidence bands and B`<sub>`MSY`</sub>` reference
- **Standard diagnostic panel** (`plot_model_fit()`): biomass, CPUE fit, residuals, harvest rate

## Installation

```r
# Development version from GitHub
# install.packages("devtools")
devtools::install_github("alistairdunn1/SurplusProductionModel")
```

### Optional dependencies for Bayesian analysis:

```r
install.packages(c("tmbstan", "rstan"))
```

## Quick Start

```r
library(SurplusProductionModel)

# Load example Ross Sea data
data(ross_sea_cpue)
data(ross_sea_catch)

# Prepare data
data_list <- list(
  cpue_data = ross_sea_cpue,
  catch_data = ross_sea_catch
)

# Fit Pella-Tomlinson model
model_fit <- fit_pella_tomlinson_model(data_list)

# Explicit state-space fit with process error
model_fit_ss <- fit_pella_tomlinson_model(
  data_list,
  options = list(process_noise = TRUE)
)

# Optional environmental covariates + AR1 process error
# (example covariate shown as annual SST anomaly)
env_data <- data.frame(
  year = ross_sea_catch$year,
  sst_anomaly = scale(rnorm(nrow(ross_sea_catch)))
)

model_fit_env_ar1 <- fit_pella_tomlinson_model(
  data = list(
    cpue_data = ross_sea_cpue,
    catch_data = ross_sea_catch,
    env_data = env_data
  ),
  options = list(
    process_noise = TRUE,
    process_error_structure = "ar1",
    env_covariates = "sst_anomaly",
    env_lag = 0
  )
)

# View results
print(model_fit)
summary(model_fit)

# Calculate reference points
ref_points <- calculate_reference_points(model_fit)
print(ref_points)

# Generate diagnostic plots
plot_model_fit(model_fit)
```

## Complete Analysis Workflow

This comprehensive example demonstrates the full analysis workflow: model fitting, convergence diagnostics, uncertainty quantification, and Bayesian inference.

### Step 1: Data Preparation and Model Fitting

```r
library(SurplusProductionModel)

# Generate synthetic example data
set.seed(123)
years <- 2010:2025
n_years <- length(years)
true_r <- 0.3
true_K <- 5000
true_m <- 2.0  # Schaefer model
true_q <- 0.001
true_B0 <- 4000
catch <- rep(400, n_years)
biomass <- numeric(n_years)
biomass[1] <- true_B0
for (t in 1:(n_years - 1)) {
  production <- true_r * biomass[t] * (1 - biomass[t] / true_K)
  biomass[t + 1] <- biomass[t] + production - catch[t]
}
cpue <- pmax(0.001, true_q * biomass * exp(rnorm(n_years, 0, 0.1)))

# Prepare data for model fitting
model_data <- list(
  cpue_data = data.frame(year = years, cpue = cpue),
  catch_data = data.frame(year = years, catch = catch)
)

# Fit the model with multi-start optimization for robustness
model_fit <- fit_pella_tomlinson_model(
  model_data, 
  n_starts = 3,  # Multiple starts to find global minimum
  verbose = TRUE
)

# View results
print(model_fit)
summary(model_fit)  # Includes parameter SEs and confidence intervals
```

### Step 2: Reference Points

```r
# Calculate biological reference points
ref_points <- calculate_reference_points(model_fit)
print(ref_points)
# Returns: MSY, BMSY, FMSY, and current B/BMSY ratio

# Calculate user-defined depletion targets analytically
target_ref_points <- calculate_reference_points(
  model_fit,
  biomass_target = 0.4,
  baseline = "K"
)
target_ref_points$target_reference_points
# Returns B_40%K and F_40%K; use baseline = "B0" when fitted B0 is available
```

### Step 3: Convergence Diagnostics

```r
# Jitter test: verify optimization found global minimum
jitter_results <- jitter_test(model_fit, n_jitter = 20, jitter_sd = 0.2)
print(jitter_results)  # Shows convergence rate and parameter stability
plot(jitter_results)   # Visualize parameter distributions across starts

# Retrospective analysis: check for systematic bias
retro_results <- retrospective_analysis(model_fit, n_peels = 5)
print(retro_results)   # Shows Mohn's rho for each parameter
plot(retro_results)    # Biomass trajectories with retrospective pattern
```

### Step 4: Model Diagnostics

```r
# Standard 4-panel diagnostic plot
plot_model_fit(model_fit)

# Enhanced residual diagnostics
plot_residuals(model_fit)  # QQ plot, histogram, residuals vs fitted

# Biomass trajectory with confidence bands
plot_biomass(model_fit)
```

### Step 5: Profile Likelihood Confidence Intervals

```r
# Profile likelihood CI for intrinsic growth rate (r)
profile_r <- profile_likelihood(model_fit, parameters = "r")
print(profile_r)
plot(profile_r)

# Profile likelihood CI for carrying capacity (K)
profile_K <- profile_likelihood(model_fit, parameters = "K")
print(profile_K)
plot(profile_K)

# Profile likelihood CI for derived quantity (MSY)
profile_msy <- profile_likelihood(
  model_fit,
  parameters = "MSY",
  n_points = 25
)

# Compact depletion-target profile example (returns B_40%K and F_40%K CIs)
profile_targets <- profile_likelihood(
  model_fit,
  parameters = c("B_40%K", "F_40%K"),
  biomass_target = 0.4,
  baseline = "K",
  n_points = 25
)
print(profile_msy)
print(profile_targets)
```

### Step 6: Bayesian Inference (Optional)

Requires `tmbstan` package for MCMC sampling via Stan's NUTS algorithm.

```r
# Fit Bayesian model (requires tmbstan)
if (requireNamespace("tmbstan", quietly = TRUE)) {
  bayes_fit <- bayesian_fit(
    model_fit,
    chains = 4,
    iter = 2000,
    warmup = 1000
  )
  print(bayes_fit)  # Summary of posterior distributions
  plot(bayes_fit)   # Posterior density plots

  # Compact depletion-target Bayesian example
  bayes_fit_targets <- bayesian_fit(
    model_fit,
    chains = 2,
    iter = 1000,
    warmup = 500,
    biomass_target = 0.4,
    baseline = "K"
  )
  subset(bayes_fit_targets$summary, parameter %in% c("B_40%K", "F_40%K"))
  
  # Posterior predictive check
  ppc <- posterior_predictive_check(bayes_fit, n_sim = 500)
  print(ppc)  # Bayesian p-value
  plot(ppc)   # Observed vs replicated data distributions
}
```

### Diagnostic Output Summary

| Function                         | Output                                                           |
| -------------------------------- | ---------------------------------------------------------------- |
| `plot_model_fit()`             | 4-panel: biomass, CPUE fit, residuals, harvest rate              |
| `plot_residuals()`             | QQ plot, histogram, residuals vs fitted                          |
| `plot_biomass()`               | Biomass trajectory with CI and B`<sub>`MSY`</sub>` reference |
| `jitter_test()`                | Optimization reliability assessment                              |
| `retrospective_analysis()`     | Mohn's rho and retrospective bias patterns                       |
| `profile_likelihood()`         | Likelihood-based confidence intervals                            |
| `bayesian_fit()`               | Full posterior distributions via MCMC                            |
| `posterior_predictive_check()` | Bayesian model validation                                        |

## Mathematical Formulation

The Pella-Tomlinson model is defined by:

**Production function:**
P(B) = r × B × (1 - (B/K)^(m-1)) / m

**State equation (deterministic default):**
B[t+1] = B[t] + P(B[t]) - C[t]

**State equation (optional state-space mode):**
log(B[t+1]) = log(B[t] + P(B[t]) - C[t]) + ε[t]

**Observation equation:**
CPUE[t] = q × B[t] × exp(η[t])

Where:

- B = biomass
- r = intrinsic growth rate
- K = carrying capacity
- m = shape parameter
- C = catch
- q = catchability coefficient
- ε ~ N(0, σ²_process) = process error in optional state-space mode (`options$process_noise = TRUE`)
- η ~ N(0, σ²_obs) = observation error

## Citation

If you use this package in your research, please cite:

```
Dunn, A. (2025). SurplusProductionModel: Pella-Tomlinson Surplus Production Model. R package version 0.1.0.
https://github.com/alistairdunn1/SurplusProductionModel
```
