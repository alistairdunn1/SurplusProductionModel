## Script to generate simulated Ross Sea toothfish dataset
## This creates a realistic dataset for demonstrating the SurplusProductionModel package
##
## Based on approximate life-history and fishery characteristics of
## Antarctic toothfish (Dissostichus mawsoni) in the Ross Sea region (CCAMLR Area 88.1/88.2)

set.seed(42)

# ===========================================================================
# True population parameters (for simulation)
# ===========================================================================
true_params <- list(
  r = 0.08, # Intrinsic growth rate (low for long-lived species)
  K = 80000, # Carrying capacity (tonnes, approximate for Ross Sea)
  m = 2.0, # Shape parameter (Schaefer model as default)
  q_LL = 0.00035, # Catchability - longline
  q_Trawl = 0.0008, # Catchability - trawl survey
  sigma_proc = 0.01, # Process error SD (low for demo identifiability)
  sigma_obs_LL = 0.15, # Observation error SD for longline CPUE
  sigma_obs_Trawl = 0.20, # Observation error SD for trawl survey CPUE
  B_initial = 72000 # Initial biomass (90% of K, pre-exploitation)
)

# ===========================================================================
# Simulation settings
# ===========================================================================
years <- 1997:2024
n_years <- length(years)

# Historical catch trajectory (tonnes)
# Fishery started in late 1990s, ramped up, was constrained, then reduced
# Note: the catch reduction in 2015-2024 creates a contrast (recovery signal)
# that helps identify production parameters — a key requirement for SPMs.
catch_trajectory <- c(
  # 1997-2002: Early development phase
  200, 500, 1100, 1500, 1800, 2200,
  # 2003-2010: Expansion and increasing catch limits
  2500, 2800, 3000, 3100, 3200, 3100, 3000, 2900,
  # 2011-2014: Peak exploitation
  2800, 2700, 2500, 2200,
  # 2015-2018: Catch limit reductions (precautionary response)
  1500, 1200, 1000, 800,
  # 2019-2024: Low catch allows partial recovery
  700, 600, 550, 500, 500, 500
)

# ===========================================================================
# Simulate population dynamics
# ===========================================================================
biomass <- numeric(n_years)
biomass[1] <- true_params$B_initial

for (t in 1:(n_years - 1)) {
  # Pella-Tomlinson production
  production <- true_params$r * biomass[t] *
    (1 - (biomass[t] / true_params$K)^(true_params$m - 1)) / true_params$m

  # Process error
  proc_error <- rnorm(1, 0, true_params$sigma_proc * biomass[t])

  # Next year biomass
  biomass[t + 1] <- max(100, biomass[t] + production - catch_trajectory[t] + proc_error)
}

# ===========================================================================
# Generate CPUE indices
# ===========================================================================

# Longline CPUE (available most years from 1998)
cpue_LL <- numeric(n_years)
cpue_LL_available <- rep(TRUE, n_years)
cpue_LL_available[1] <- FALSE # No longline CPUE in first year

for (t in 1:n_years) {
  if (cpue_LL_available[t]) {
    cpue_LL[t] <- true_params$q_LL * biomass[t] *
      exp(rnorm(1, -0.5 * true_params$sigma_obs_LL^2, true_params$sigma_obs_LL))
  } else {
    cpue_LL[t] <- NA
  }
}

# Trawl survey CPUE (biennial survey from 2000, with some gaps)
cpue_Trawl <- rep(NA_real_, n_years)
survey_years <- seq(2000, 2024, by = 2)
# Remove a couple of years to simulate survey gaps
survey_years <- survey_years[!survey_years %in% c(2020)] # No survey in COVID year

for (yr in survey_years) {
  t <- which(years == yr)
  if (length(t) == 1) {
    cpue_Trawl[t] <- true_params$q_Trawl * biomass[t] *
      exp(rnorm(1, -0.5 * true_params$sigma_obs_Trawl^2, true_params$sigma_obs_Trawl))
  }
}

# ===========================================================================
# Assemble data frames following package schema
# ===========================================================================

# Catch data (single area)
ross_sea_catch <- data.frame(
  year = years,
  catch = catch_trajectory,
  stringsAsFactors = FALSE
)
attr(ross_sea_catch, "units") <- "tonnes"
attr(ross_sea_catch, "source") <- "Simulated (based on CCAMLR Area 88.1/88.2 characteristics)"
attr(ross_sea_catch, "created") <- Sys.time()

# CPUE data — longline index
cpue_ll_df <- data.frame(
  year = years[cpue_LL_available],
  cpue = cpue_LL[cpue_LL_available],
  label = "Longline",
  stringsAsFactors = FALSE
)

# CPUE data — trawl survey index
trawl_idx <- !is.na(cpue_Trawl)
cpue_trawl_df <- data.frame(
  year = years[trawl_idx],
  cpue = cpue_Trawl[trawl_idx],
  label = "Trawl",
  stringsAsFactors = FALSE
)

# Combined CPUE data (multi-index)
ross_sea_cpue <- rbind(cpue_ll_df, cpue_trawl_df)
ross_sea_cpue <- ross_sea_cpue[order(ross_sea_cpue$year, ross_sea_cpue$label), ]
rownames(ross_sea_cpue) <- NULL
attr(ross_sea_cpue, "units") <- "kg/hook (Longline), kg/km2 (Trawl)"
attr(ross_sea_cpue, "source") <- "Simulated (based on CCAMLR Area 88.1/88.2 characteristics)"
attr(ross_sea_cpue, "created") <- Sys.time()

# Also provide a simple single-index CPUE (longline only) for basic examples
ross_sea_cpue_simple <- cpue_ll_df[, c("year", "cpue")]
rownames(ross_sea_cpue_simple) <- NULL
attr(ross_sea_cpue_simple, "units") <- "kg/hook"
attr(ross_sea_cpue_simple, "source") <- "Simulated (based on CCAMLR Area 88.1/88.2 characteristics)"
attr(ross_sea_cpue_simple, "created") <- Sys.time()

# True parameter values (for validation)
ross_sea_true_params <- true_params

# True biomass trajectory (for validation)
ross_sea_true_biomass <- data.frame(
  year = years,
  biomass = biomass
)

# ===========================================================================
# Save datasets
# ===========================================================================
usethis::use_data(ross_sea_catch, overwrite = TRUE)
usethis::use_data(ross_sea_cpue, overwrite = TRUE)
usethis::use_data(ross_sea_cpue_simple, overwrite = TRUE)
usethis::use_data(ross_sea_true_params, overwrite = TRUE)
usethis::use_data(ross_sea_true_biomass, overwrite = TRUE)

cat("Datasets created successfully:\n")
cat("  ross_sea_catch         - Catch data (1997-2024)\n")
cat("  ross_sea_cpue          - Multi-index CPUE (Longline + Trawl survey)\n")
cat("  ross_sea_cpue_simple   - Single-index CPUE (Longline only)\n")
cat("  ross_sea_true_params   - True simulation parameters (for validation)\n")
cat("  ross_sea_true_biomass  - True biomass trajectory (for validation)\n")
