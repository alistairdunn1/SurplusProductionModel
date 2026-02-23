# SurplusProductionModel

[![R-CMD-check](https://github.com/alistairdunn1/SurplusProductionModel/workflows/R-CMD-check/badge.svg)](https://github.com/alistairdunn1/SurplusProductionModel/actions)
[![codecov](https://codecov.io/gh/alistairdunn1/SurplusProductionModel/branch/main/graph/badge.svg)](https://codecov.io/gh/alistairdunn1/SurplusProductionModel)
[![pkgdown](https://github.com/alistairdunn1/SurplusProductionModel/actions/workflows/pkgdown.yaml/badge.svg)](https://alistairdunn1.github.io/SurplusProductionModel/)

## Overview

`SurplusProductionModel` is an R package that implements a state-space spatial Pella-Tomlinson surplus production model for Antarctic toothfish (*Dissostichus mawsoni*) stock assessment. This package serves as the operating model foundation for the ATO rTMB project's Management Strategy Evaluation (MSE) framework.

## Key Features

### Model Fitting
- **Pella-Tomlinson production function** with flexible shape parameter (*m*)
  - Schaefer model (*m* = 2): symmetric production curve
  - Fox model (*m* = 1): asymmetric production curve with peak at lower biomass
- **State-space framework** with separate process and observation error
- **RTMB integration** for automatic differentiation and fast optimization via `nlminb`
- **Spatial structure** supporting multiple management areas with optional movement
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
profile_r <- profile_likelihood(model_fit, parameter = "r")
print(profile_r)
plot(profile_r)

# Profile likelihood CI for carrying capacity (K)
profile_K <- profile_likelihood(model_fit, parameter = "K")
print(profile_K)
plot(profile_K)

# Profile likelihood CI for derived quantity (MSY)
profile_msy <- profile_likelihood(
  model_fit, 
  derived_quantity = "MSY",
  n_points = 25
)
print(profile_msy)
plot(profile_msy)
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
  
  # Posterior predictive check
  ppc <- posterior_predictive_check(bayes_fit, n_sim = 500)
  print(ppc)  # Bayesian p-value
  plot(ppc)   # Observed vs replicated data distributions
}
```

### Diagnostic Output Summary

| Function | Output |
|----------|--------|
| `plot_model_fit()` | 4-panel: biomass, CPUE fit, residuals, harvest rate |
| `plot_residuals()` | QQ plot, histogram, residuals vs fitted |
| `plot_biomass()` | Biomass trajectory with CI and B<sub>MSY</sub> reference |
| `jitter_test()` | Optimization reliability assessment |
| `retrospective_analysis()` | Mohn's rho and retrospective bias patterns |
| `profile_likelihood()` | Likelihood-based confidence intervals |
| `bayesian_fit()` | Full posterior distributions via MCMC |
| `posterior_predictive_check()` | Bayesian model validation |

## Mathematical Formulation

The Pella-Tomlinson model is defined by:

**Production function:**
P(B) = r × B × (1 - (B/K)^(m-1)) / m

**State equation:**
B[t+1] = B[t] + P(B[t]) - C[t] + ε[t]

**Observation equation:**
CPUE[t] = q × B[t] × exp(η[t])

Where:
- B = biomass
- r = intrinsic growth rate
- K = carrying capacity  
- m = shape parameter
- C = catch
- q = catchability coefficient
- ε ~ N(0, σ²_process) = process error
- η ~ N(0, σ²_obs) = observation error

## Project Context

This package is part of the larger ATO rTMB (Antarctic Toothfish Assessment with RTMB) project, which aims to develop a comprehensive suite of stock assessment and management strategy evaluation tools for Antarctic toothfish.

**Related packages:**
- MSE: Management Strategy Evaluation framework (depends on this package)
- IntegratedAgelengthModel: Comprehensive age-length structured assessment (future development)

## Contributing

Please read our [contributing guidelines](CONTRIBUTING.md) and [code of conduct](CODE_OF_CONDUCT.md).

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Citation

If you use this package in your research, please cite:

```
Dunn, A. (2025). SurplusProductionModel: Pella-Tomlinson Surplus Production Model 
for Antarctic Toothfish. R package version 0.1.0.
https://github.com/alistairdunn1/SurplusProductionModel
```

## Support

- Report bugs: [GitHub Issues](https://github.com/alistairdunn1/SurplusProductionModel/issues)
- Documentation: [Package website](https://alistairdunn1.github.io/SurplusProductionModel/)
- Questions: [Discussions](https://github.com/alistairdunn1/SurplusProductionModel/discussions)
