#' SurplusProductionModel: Pella-Tomlinson Surplus Production Model for Antarctic Toothfish
#'
#' This package implements a state-space spatial Pella-Tomlinson surplus production 
#' model for Antarctic toothfish (Dissostichus mawsoni) stock assessment. It serves 
#' as the MVP foundation for the ATO rTMB project.
#'
#' @section Key Features:
#' 
#' \itemize{
#'   \item Pella-Tomlinson production function with flexible shape parameter (m)
#'   \item Special cases: Schaefer model (m=2) and Fox model (m=1)
#'   \item State-space framework with process and observation error
#'   \item RTMB integration for automatic differentiation and optimization
#'   \item Spatial structure supporting multiple management areas
#'   \item Reference point calculations (MSY, BMSY, FMSY)
#'   \item Comprehensive diagnostics and model validation tools
#' }
#'
#' @section Main Functions:
#' 
#' \itemize{
#'   \item \code{\link{fit_pella_tomlinson_model}}: Fit the surplus production model
#'   \item \code{\link{calculate_reference_points}}: Calculate biological reference points
#'   \item \code{\link{estimate_biomass}}: Extract biomass estimates
#'   \item \code{\link{plot_model_fit}}: Generate diagnostic plots
#' }
#'
#' @section Mathematical Formulation:
#' 
#' The Pella-Tomlinson model is defined by:
#' 
#' Production function: P(B) = r × B × (1 - (B/K)^(m-1)) / m
#' 
#' State equation: B[t+1] = B[t] + P(B[t]) - C[t] + ε[t]
#' 
#' Observation equation: CPUE[t] = q × B[t] × exp(η[t])
#'
#' @docType package
#' @name SurplusProductionModel-package
#' @aliases SurplusProductionModel
#' @keywords package
NULL