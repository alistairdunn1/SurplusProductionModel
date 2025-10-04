#' Diagnostic Plot Functions for Surplus Production Model
#'
#' This file contains diagnostic plotting functions for the Pella-Tomlinson surplus production model.
#'
#' @name diagnostics
NULL

#' Plot Model Fit Diagnostics
#'
#' Generates diagnostic plots for a fitted ProductionModel object:
#' - Biomass trajectory
#' - CPUE observed vs fitted
#' - Residuals
#' - Harvest rate
#'
#' @param model_fit A fitted ProductionModel object
#' @param which Character vector specifying which plots to show ("biomass", "cpue", "residuals", "harvest", "all")
#' @param show Logical, whether to print plots (default: TRUE)
#' @return List of ggplot objects
#' @export
plot_model_fit <- function(model_fit, which = "all", show = TRUE) {
  if (!inherits(model_fit, "ProductionModel")) stop("Input must be a ProductionModel object")
  if (!model_fit@fitted) stop("Model must be fitted before plotting")
  results <- model_fit@results
  years <- model_fit@data$years
  plots <- list()

  # Biomass trajectory
  if ("biomass" %in% which || which == "all") {
    df <- data.frame(year = years, biomass = results$biomass)
    p_biomass <- ggplot2::ggplot(df, ggplot2::aes(x = year, y = biomass)) +
      ggplot2::geom_line(color = "#0072B2", size = 1.2) +
      ggplot2::geom_point(color = "#0072B2") +
      ggplot2::labs(title = "Estimated Biomass Trajectory", y = "Biomass (tonnes)", x = "Year") +
      ggplot2::theme_bw()
    plots$biomass <- p_biomass
    if (show) print(p_biomass)
  }

  # CPUE observed vs fitted
  if (("cpue" %in% which || which == "all") && "fitted_cpue" %in% names(results)) {
    df <- data.frame(year = years, observed = model_fit@data$cpue, fitted = results$fitted_cpue)
    p_cpue <- ggplot2::ggplot(df) +
      ggplot2::geom_point(ggplot2::aes(x = year, y = observed), color = "#D55E00", size = 2) +
      ggplot2::geom_line(ggplot2::aes(x = year, y = fitted), color = "#009E73", size = 1.2) +
      ggplot2::labs(title = "CPUE: Observed vs Fitted", y = "CPUE", x = "Year") +
      ggplot2::theme_bw()
    plots$cpue <- p_cpue
    if (show) print(p_cpue)
  }

  # Residuals plot
  if (("residuals" %in% which || which == "all") && "residuals" %in% names(results)) {
    df <- data.frame(year = years, residuals = results$residuals)
    p_resid <- ggplot2::ggplot(df, ggplot2::aes(x = year, y = residuals)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
      ggplot2::geom_line(color = "#CC79A7", size = 1.2) +
      ggplot2::geom_point(color = "#CC79A7") +
      ggplot2::labs(title = "Log-Residuals: Observed - Fitted CPUE", y = "Log-Residual", x = "Year") +
      ggplot2::theme_bw()
    plots$residuals <- p_resid
    if (show) print(p_resid)
  }

  # Harvest rate plot
  if (("harvest" %in% which || which == "all") && "harvest_rate" %in% names(results)) {
    df <- data.frame(year = years, harvest_rate = results$harvest_rate)
    p_harvest <- ggplot2::ggplot(df, ggplot2::aes(x = year, y = harvest_rate)) +
      ggplot2::geom_line(color = "#F0E442", size = 1.2) +
      ggplot2::geom_point(color = "#F0E442") +
      ggplot2::labs(title = "Harvest Rate Trajectory", y = "Harvest Rate (Catch/Biomass)", x = "Year") +
      ggplot2::theme_bw()
    plots$harvest <- p_harvest
    if (show) print(p_harvest)
  }

  invisible(plots)
}
