#' SurplusProductionModel: Pella-Tomlinson Surplus Production Model
#'
#' This package implements a spatial Pella-Tomlinson surplus production model. It also supports
#' an optional state-space formulation with process error for model fitting. It
#' serves as the operating model foundation for the ATO rTMB project's
#' Management Strategy Evaluation (MSE) framework.
#'
#' @section Key Features:
#'
#' \itemize{
#'   \item Pella-Tomlinson production function with flexible shape parameter (m)
#'   \item Special cases: Schaefer model (m=2) and Fox model (m=1)
#'   \item Optional state-space framework with process error, plus observation error in CPUE
#'   \item RTMB integration for automatic differentiation and optimization
#'   \item Spatial structure supporting multiple management areas with movement
#'   \item Multi-index support for multiple CPUE series per area
#'   \item Convergence diagnostics: multi-start optimization, jitter tests, retrospective analysis
#'   \item Profile likelihood confidence intervals for parameters and derived quantities
#'   \item Bayesian inference via tmbstan MCMC sampling
#'   \item Comprehensive diagnostics and model validation tools
#' }
#'
#' @section Model Fitting:
#'
#' \itemize{
#'   \item \code{\link{fit_pella_tomlinson_model}}: Fit the surplus production model
#'   \item \code{\link{calculate_reference_points}}: Calculate biological reference points (MSY, BMSY, FMSY)
#'   \item \code{\link{estimate_biomass}}: Extract biomass estimates with confidence intervals
#' }
#'
#' @section Convergence Diagnostics:
#'
#' \itemize{
#'   \item \code{\link{jitter_test}}: Test optimization reliability from perturbed starting values
#'   \item \code{\link{retrospective_analysis}}: Compute Mohn's rho for retrospective bias assessment
#' }
#'
#' @section Uncertainty Quantification:
#'
#' \itemize{
#'   \item \code{\link{profile_likelihood}}: Profile likelihood CIs for parameters and derived quantities
#'   \item \code{\link{bayesian_fit}}: Bayesian MCMC sampling via tmbstan
#'   \item \code{\link{posterior_predictive_check}}: Bayesian model validation
#' }
#'
#' @section Diagnostics and Visualization:
#'
#' \itemize{
#'   \item \code{\link{plot_model_fit}}: Standard 4-panel diagnostic plots
#'   \item \code{\link{plot_residuals}}: Enhanced residual diagnostics (QQ, histogram, fitted)
#'   \item \code{\link{plot_biomass}}: Biomass trajectories with confidence bands
#' }
#'
#' @section Mathematical Formulation:
#'
#' The Pella-Tomlinson model is defined by:
#'
#' Production function: P(B) = r × B × (1 - (B/K)^(m-1)) / m
#'
#' State equation: \code{B(t+1) = B(t) + P(B(t)) - C(t) + e(t)}
#'
#' Observation equation: \code{CPUE(t) = q * B(t) * exp(n(t))}
#'
"_PACKAGE"
