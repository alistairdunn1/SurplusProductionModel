# SurplusProductionModel

[![R-CMD-check](https://github.com/alistairdunn1/SurplusProductionModel/workflows/R-CMD-check/badge.svg)](https://github.com/alistairdunn1/SurplusProductionModel/actions)
[![codecov](https://codecov.io/gh/alistairdunn1/SurplusProductionModel/branch/main/graph/badge.svg)](https://codecov.io/gh/alistairdunn1/SurplusProductionModel)

## Overview

`SurplusProductionModel` is an R package that implements a state-space spatial Pella-Tomlinson surplus production model for Antarctic toothfish (*Dissostichus mawsoni*) stock assessment. This package serves as the MVP foundation for the ATO rTMB project.

## Key Features

- **Pella-Tomlinson production function** with flexible shape parameter (m)
- **Special cases**: Schaefer model (m=2) and Fox model (m=1)
- **State-space framework** with process and observation error
- **RTMB integration** for automatic differentiation and optimization
- **Spatial structure** supporting multiple management areas
- **Reference point calculations** (MSY, BMSY, FMSY)
- **Comprehensive diagnostics** and model validation tools

## Installation

```r
# Development version from GitHub
# install.packages("devtools")
devtools::install_github("alistairdunn1/SurplusProductionModel")
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

# Calculate reference points
ref_points <- calculate_reference_points(model_fit)
print(ref_points)

# Generate diagnostic plots
plot_model_fit(model_fit)
```

## Worked Example: Model Fitting and Diagnostics

This example demonstrates fitting the Pella-Tomlinson surplus production model to synthetic data and generating diagnostic plots.

```r
library(SurplusProductionModel)

# Generate synthetic example data
set.seed(123)
years <- 2010:2020
n_years <- length(years)
true_r <- 0.3
true_K <- 5000
true_m <- 2.0  # Schaefer model
true_q <- 0.001
true_B0 <- 4000
catch <- rep(800, n_years)
biomass <- numeric(n_years)
biomass[1] <- true_B0
for (t in 1:(n_years - 1)) {
  production <- true_r * biomass[t] * (1 - biomass[t] / true_K)
  biomass[t + 1] <- biomass[t] + production - catch[t]
}
cpue <- pmax(0.001, true_q * biomass * exp(rnorm(n_years, 0, 0.05)))

# Prepare data for model fitting
model_data <- list(
  cpue_data = data.frame(year = years, cpue = cpue),
  catch_data = data.frame(year = years, catch = catch)
)

# Fit the model
model_fit <- fit_pella_tomlinson_model(model_data)

# Print fitted parameters
print(model_fit)

# Calculate reference points
ref_points <- calculate_reference_points(model_fit)
print(ref_points)

# Generate diagnostic plots
plot_model_fit(model_fit)
```

**Diagnostic plots produced:**
- Estimated biomass trajectory
- CPUE: observed vs fitted
- Log-residuals (CPUE)
- Harvest rate trajectory

These plots help assess model fit, residual structure, and exploitation trends.

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
