# SurplusProductionModel

## Overview

`SurplusProductionModel` is an R package that implements a spatial Pella-Tomlinson surplus production model. It supports an optional state-space formulation with process error for model fitting. This package also serves as the operating model foundation for the Surplus Production Model Management Strategy Evaluation (MSE) framework.

## Key Features

### Model Fitting

- **Pella-Tomlinson production function** with flexible shape parameter (*m*) (Pella and Tomlinson 1969)
  - Schaefer model (*m* = 2): symmetric production curve (Schaefer 1954)
  - Fox model (*m* = 1): asymmetric production curve with peak at lower biomass (Fox 1970)
- **Optional state-space framework** with process error in biomass dynamics when `options$process_noise = TRUE`, plus observation error in CPUE
- **RTMB integration** for automatic differentiation and fast optimization via `nlminb`
- **Spatial structure** with independently estimated carrying capacities per area and optional movement via a gravity kernel or a complete transition matrix
- **Multi-index support** for multiple CPUE series per area, or a shared catchability across areas for unlabelled indices

### Convergence Diagnostics

- **Multi-start optimization** (`options$n_starts`) to search for a better optimum
- **Jitter test** (`jitter_test()`) to assess optimization reliability from perturbed starts
- **Retrospective analysis** (`retrospective_analysis()`) with Mohn's rho for systematic bias detection

### Uncertainty Quantification

- **Delta-method standard errors** via `sdreport` (Hessian-based)
- **Profile likelihood confidence intervals** (`profile_likelihood()`) for parameters and derived quantities
- **Bayesian MCMC sampling** (`bayesian_fit()`) via the sparse No-U-Turn Sampler in `SparseNUTS`
- **Posterior predictive checks** (`posterior_predictive_check()`) for model validation

### Diagnostics & Visualization

- **Reference points**: MSY, B<sub>MSY</sub>, F<sub>MSY</sub> with current B/B<sub>MSY</sub> ratio
- **Residual diagnostics** (`plot_residuals()`): QQ plots, residuals vs. fitted, histograms
- **Biomass plots** (`plot_biomass()`): trajectories with confidence bands and B<sub>MSY</sub> reference
- **Standard diagnostic panel** (`plot_model_fit()`): biomass, CPUE fit, residuals, harvest rate

## Installation

```r
# Development version from GitHub
# install.packages("devtools")
devtools::install_github("alistairdunn1/SurplusProductionModel")
```

### Optional dependency for Bayesian analysis:

```r
# install.packages("remotes")
remotes::install_github("noaa-afsc/SparseNUTS")
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

# Fit the Schaefer special case (fix the poorly identified shape parameter)
model_fit <- fit_pella_tomlinson_model(
  data_list, options = list(fixed_params = list(log_m = log(2)))
)

# Explicit state-space fit with process error
model_fit_ss <- fit_pella_tomlinson_model(
  data_list,
  options = list(process_noise = TRUE, fixed_params = list(log_m = log(2)))
)

# Optional environmental covariates + AR1 process error
# (example covariate shown as annual SST anomaly)
set.seed(123)
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
    fixed_params = list(log_m = log(2)),
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

## Spatial Models and Fixed Fox Fits

Add an `area` column to catch and CPUE data for spatial fits. Carrying capacity
is estimated independently per area (`area_k_shares` is deprecated). By default,
catchability is estimated per area, or per area and label if CPUE has a `label`
column. Set `options$shared_q = TRUE` for one catchability across two or more
areas: omit the `label` column and provide explicit `params_init` with
`log_q_shared` replacing all per-area `log_q.<area>` entries. Automatic starts
are not supported in this mode.

Movement accepts either gravity inputs (`distance_matrix`, `attractiveness`,
`decay`, and `movement_rate`) or a complete `data$movement$transition_matrix`.
The latter has source areas in rows and destinations in columns, with finite,
non-negative entries and each row summing to one, including retention. It is
applied in full and overrides gravity inputs; `movement_rate` is ignored and
cannot be estimated with a complete matrix. Movement inputs trigger a 50-year
zero-catch spin-up when `spinup_years` is zero. Initial depletion multiplies
the resulting unfished spatial equilibrium, which can differ from per-area `K`.

For the exact Fox model, use `options = list(fixed_params = list(log_m = 0))`.
The objective then uses the analytic limit at *m* = 1.

See the [getting-started vignette](vignettes/getting-started.Rmd),
[multi-index example](inst/examples/example_multiindex.R), and
[shared catchability and movement example](inst/examples/example_shared_q.R).

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
true_B_initial <- 4000
catch <- rep(400, n_years)
biomass <- numeric(n_years)
biomass[1] <- true_B_initial
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
  options = list(
    n_starts = 3,  # Compare multiple starting points
    control = list(eval.max = 5000, iter.max = 2000)
  )
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
# Returns B_40%K and F_40%K; use baseline = "B_initial" for B_initial = d0 * K
```

### Step 3: Convergence Diagnostics

```r
# Jitter test: assess sensitivity to starting values
jitter_results <- jitter_test(model_fit, n_jitter = 20, jitter_sd = 0.2)
print(jitter_results)  # Shows convergence rate and parameter stability
plot(jitter_results)   # Visualize parameter distributions across starts

# Retrospective analysis: check for systematic bias
retro_results <- retrospective_analysis(model_fit, n_peels = 5)
print(retro_results)   # Shows Mohn's rho for terminal biomass
plot(retro_results)    # Biomass trajectories with retrospective pattern
```

### Step 4: Model Diagnostics

Known plotting limitation: single-index CPUE observations and fitted values
can appear in separate `obs_mat` and `A1` facets despite representing the same
area. This is a plotting issue.

```r
# Standard 4-panel diagnostic plot
plot_model_fit(model_fit)

# Enhanced residual diagnostics
plot_residuals(model_fit)  # QQ plot, histogram, residuals vs fitted

# Biomass trajectory with confidence bands
plot_biomass(model_fit)
```

### Step 5: Profile Likelihood Confidence Intervals

The workflow above estimates *m*. Derived-quantity profiles currently require
*m* to be estimated: fixed-shape fits fail when profiling MSY or depletion
targets. An `NA` confidence limit means the evaluated profile grid did not
establish that limit.

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

Requires the `SparseNUTS` package for MCMC sampling via its sparse NUTS
implementation.

```r
# Fit Bayesian model (requires SparseNUTS)
if (requireNamespace("SparseNUTS", quietly = TRUE)) {
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
  ppc <- posterior_predictive_check(bayes_fit, n_sims = 500)
  print(ppc)  # Bayesian p-value
  plot(ppc)   # Observed vs replicated data distributions
}
```

### Diagnostic Output Summary

| Function                         | Output                                                             |
| -------------------------------- | ------------------------------------------------------------------ |
| `plot_model_fit()`             | 4-panel: biomass, CPUE fit, residuals, harvest rate                |
| `plot_residuals()`             | QQ plot, histogram, residuals vs fitted                            |
| `plot_biomass()`               | Biomass trajectory with CI and B<sub>MSY</sub> reference |
| `jitter_test()`                | Optimization reliability assessment                                |
| `retrospective_analysis()`     | Mohn's rho and retrospective bias patterns                         |
| `profile_likelihood()`         | Likelihood-based confidence intervals                              |
| `bayesian_fit()`               | Full posterior distributions via MCMC                              |
| `posterior_predictive_check()` | Bayesian model validation                                          |

## Mathematical Formulation

The Pella-Tomlinson model is defined by:

**Production function (standard Pella-Tomlinson):**
P(B) = r / (m - 1) × B × (1 - (B/K)^(m-1))

with the Fox limit P(B) = r × B × log(K/B) as *m* → 1, and the Schaefer
form P(B) = r × B × (1 - B/K) at *m* = 2. Here *r* scales productivity;
the low-biomass per-capita limit is r/(m-1) for m > 1 and is unbounded in
the Fox limit. Reference points are B_MSY = K·*m*^(-1/(m-1)), F_MSY = *r*/*m*,
and MSY = F_MSY·B_MSY (Fox: K/e, *r*, rK/e).

**State equation (deterministic default, without movement or covariates):**
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
Dunn, A. (2026). SurplusProductionModel: Pella-Tomlinson Surplus Production Model. R package version 0.1.3.
https://github.com/alistairdunn1/SurplusProductionModel
```
