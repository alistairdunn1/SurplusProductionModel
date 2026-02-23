#' Diagnostic Plot Functions for Surplus Production Model
#'
#' This file contains diagnostic plotting functions for the Pella-Tomlinson surplus production model.
#'
#' @name diagnostics
NULL

# Internal helper: melt a named array to a long data.frame (replaces reshape2::melt)
melt_array <- function(arr, varnames, value_name) {
  dn <- dimnames(arr)
  grid <- do.call(expand.grid, c(lapply(dn, function(x) x), list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)))
  names(grid) <- varnames
  grid[[value_name]] <- as.vector(arr)
  grid
}

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
  if (!model_fit$fitted) stop("Model must be fitted before plotting")
  results <- model_fit$results
  years <- model_fit$data$years
  plots <- list()

  # Biomass trajectory
  if ("biomass" %in% which || which == "all") {
    biomass <- results$biomass
    if (is.matrix(biomass)) {
      df <- as.data.frame(biomass)
      df$year <- years
      df_long <- tidyr::pivot_longer(df, cols = -year, names_to = "area", values_to = "biomass")
      p_biomass <- ggplot2::ggplot(df_long, ggplot2::aes(x = year, y = biomass)) +
        ggplot2::geom_line(color = "#0072B2", size = 1.0) +
        ggplot2::geom_point(color = "#0072B2", size = 1.5) +
        ggplot2::facet_wrap(~area, scales = "free_y") +
        ggplot2::labs(title = "Estimated Biomass Trajectory", y = "Biomass (tonnes)", x = "Year") +
        ggplot2::theme_bw()
    } else {
      df <- data.frame(year = years, biomass = biomass)
      p_biomass <- ggplot2::ggplot(df, ggplot2::aes(x = year, y = biomass)) +
        ggplot2::geom_line(color = "#0072B2", size = 1.2) +
        ggplot2::geom_point(color = "#0072B2") +
        ggplot2::labs(title = "Estimated Biomass Trajectory", y = "Biomass (tonnes)", x = "Year") +
        ggplot2::theme_bw()
    }
    plots$biomass <- p_biomass
    if (show) print(p_biomass)
  }

  # CPUE observed vs fitted
  if (("cpue" %in% which || which == "all") && "fitted_cpue" %in% names(results)) {
    obs <- model_fit$data$cpue
    fit <- results$fitted_cpue
    # If multiple indices per area, plot per label (facet by area and label)
    if (is.array(obs) && length(dim(obs)) == 3) {
      # build long data frame for observed and fitted
      dn <- dimnames(obs)
      label_vals <- dn[[3]]
      # Melt arrays
      obs_long <- melt_array(obs, c("year", "area", "label"), "observed")
      if (is.array(fit) && length(dim(fit)) == 3) {
        fit_long <- melt_array(fit, c("year", "area", "label"), "fitted")
      } else {
        # if fitted is matrix [year x area] but obs has labels, recycle across labels
        fit_mat <- if (is.matrix(fit)) fit else matrix(fit, ncol = 1)
        fit_arr <- array(rep(fit_mat, length(label_vals)), dim = c(nrow(fit_mat), ncol(fit_mat), length(label_vals)), dimnames = list(year = rownames(fit_mat), area = colnames(fit_mat), label = label_vals))
        fit_long <- melt_array(fit_arr, c("year", "area", "label"), "fitted")
      }
      df_long <- merge(obs_long, fit_long, by = c("year", "area", "label"), all = TRUE)
      # coerce year to integer
      df_long$year <- as.integer(as.character(df_long$year))
      p_cpue <- ggplot2::ggplot(df_long) +
        ggplot2::geom_point(ggplot2::aes(x = year, y = observed), color = "#D55E00", size = 1.2, alpha = 0.7) +
        ggplot2::geom_line(ggplot2::aes(x = year, y = fitted), color = "#009E73", size = 0.9) +
        ggplot2::facet_grid(label ~ area, scales = "free_y") +
        ggplot2::labs(title = "CPUE: Observed vs Fitted (per index)", y = "CPUE", x = "Year") +
        ggplot2::theme_bw()
    } else {
      obs_mat <- obs
      if (is.matrix(obs_mat) || is.matrix(fit)) {
        fit_mat <- if (is.matrix(fit)) fit else matrix(fit, ncol = 1)
        df <- as.data.frame(obs_mat)
        df$year <- years
        df_obs <- tidyr::pivot_longer(df, cols = -year, names_to = "area", values_to = "observed")
        df2 <- as.data.frame(fit_mat)
        df2$year <- years
        df_fit <- tidyr::pivot_longer(df2, cols = -year, names_to = "area", values_to = "fitted")
        df_long <- merge(df_obs, df_fit, by = c("year", "area"), all = TRUE)
        p_cpue <- ggplot2::ggplot(df_long) +
          ggplot2::geom_point(ggplot2::aes(x = year, y = observed), color = "#D55E00", size = 1.5, alpha = 0.8) +
          ggplot2::geom_line(ggplot2::aes(x = year, y = fitted), color = "#009E73", size = 1.0) +
          ggplot2::facet_wrap(~area, scales = "free_y") +
          ggplot2::labs(title = "CPUE: Observed vs Fitted", y = "CPUE", x = "Year") +
          ggplot2::theme_bw()
      } else {
        df <- data.frame(year = years, observed = obs, fitted = fit)
        p_cpue <- ggplot2::ggplot(df) +
          ggplot2::geom_point(ggplot2::aes(x = year, y = observed), color = "#D55E00", size = 2) +
          ggplot2::geom_line(ggplot2::aes(x = year, y = fitted), color = "#009E73", size = 1.2) +
          ggplot2::labs(title = "CPUE: Observed vs Fitted", y = "CPUE", x = "Year") +
          ggplot2::theme_bw()
      }
    }
    plots$cpue <- p_cpue
    if (show) print(p_cpue)
  }

  # Residuals plot
  if (("residuals" %in% which || which == "all") && "residuals" %in% names(results)) {
    res <- results$residuals
    if (is.matrix(res)) {
      df <- as.data.frame(res)
      df$year <- years
      df_long <- tidyr::pivot_longer(df, cols = -year, names_to = "area", values_to = "residuals")
      p_resid <- ggplot2::ggplot(df_long, ggplot2::aes(x = year, y = residuals)) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
        ggplot2::geom_line(color = "#CC79A7", size = 1.0) +
        ggplot2::geom_point(color = "#CC79A7", size = 1.5) +
        ggplot2::facet_wrap(~area) +
        ggplot2::labs(title = "Log-Residuals: Observed - Fitted CPUE", y = "Log-Residual", x = "Year") +
        ggplot2::theme_bw()
    } else if (is.array(res) && length(dim(res)) == 3) {
      # facet by area and label when residuals are per index
      res_long <- melt_array(res, c("year", "area", "label"), "residuals")
      res_long$year <- as.integer(as.character(res_long$year))
      p_resid <- ggplot2::ggplot(res_long, ggplot2::aes(x = year, y = residuals)) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
        ggplot2::geom_line(color = "#CC79A7", size = 0.9) +
        ggplot2::geom_point(color = "#CC79A7", size = 1.2) +
        ggplot2::facet_grid(label ~ area) +
        ggplot2::labs(title = "Log-Residuals by Index", y = "Log-Residual", x = "Year") +
        ggplot2::theme_bw()
    } else {
      df <- data.frame(year = years, residuals = res)
      p_resid <- ggplot2::ggplot(df, ggplot2::aes(x = year, y = residuals)) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
        ggplot2::geom_line(color = "#CC79A7", size = 1.2) +
        ggplot2::geom_point(color = "#CC79A7") +
        ggplot2::labs(title = "Log-Residuals: Observed - Fitted CPUE", y = "Log-Residual", x = "Year") +
        ggplot2::theme_bw()
    }
    plots$residuals <- p_resid
    if (show) print(p_resid)
  }

  # Harvest rate plot
  if (("harvest" %in% which || which == "all") && "harvest_rate" %in% names(results)) {
    hr <- results$harvest_rate
    if (is.matrix(hr)) {
      df <- as.data.frame(hr)
      df$year <- years
      df_long <- tidyr::pivot_longer(df, cols = -year, names_to = "area", values_to = "harvest_rate")
      p_harvest <- ggplot2::ggplot(df_long, ggplot2::aes(x = year, y = harvest_rate)) +
        ggplot2::geom_line(color = "#F0E442", size = 1.0) +
        ggplot2::geom_point(color = "#F0E442", size = 1.5) +
        ggplot2::facet_wrap(~area) +
        ggplot2::labs(title = "Harvest Rate Trajectory", y = "Harvest Rate (Catch/Biomass)", x = "Year") +
        ggplot2::theme_bw()
    } else {
      df <- data.frame(year = years, harvest_rate = hr)
      p_harvest <- ggplot2::ggplot(df, ggplot2::aes(x = year, y = harvest_rate)) +
        ggplot2::geom_line(color = "#F0E442", size = 1.2) +
        ggplot2::geom_point(color = "#F0E442") +
        ggplot2::labs(title = "Harvest Rate Trajectory", y = "Harvest Rate (Catch/Biomass)", x = "Year") +
        ggplot2::theme_bw()
    }
    plots$harvest <- p_harvest
    if (show) print(p_harvest)
  }

  invisible(plots)
}


#' Plot Residual Diagnostics
#'
#' Generates detailed residual diagnostic plots including QQ plots,
#' residuals vs fitted, and residual histograms.  For multi-area or
#' multi-index models the plots are faceted by area and/or label.
#'
#' @param model_fit A fitted ProductionModel object.
#' @param which Character vector: one or more of \code{"qq"}, \code{"fitted"},
#'   \code{"histogram"}, or \code{"all"} (default).
#' @param show Logical, whether to print plots (default: TRUE).
#' @return List of ggplot objects (invisibly).
#' @export
plot_residuals <- function(model_fit, which = "all", show = TRUE) {
  if (!inherits(model_fit, "ProductionModel")) stop("Input must be a ProductionModel object")
  if (!model_fit$fitted) stop("Model must be fitted before diagnostic plotting")

  results <- model_fit$results
  years   <- model_fit$data$years
  plots   <- list()

  # Build a long-form residuals data.frame
  res     <- results$residuals
  fit_val <- results$fitted_cpue

  if (is.array(res) && length(dim(res)) == 3) {
    df <- melt_array(res, c("year", "area", "label"), "residual")
    fit_df <- melt_array(fit_val, c("year", "area", "label"), "fitted_val")
    df <- merge(df, fit_df, by = c("year", "area", "label"))
    df$year <- as.integer(as.character(df$year))
    facet_fn <- function(p) p + ggplot2::facet_grid(label ~ area)
  } else if (is.matrix(res)) {
    df <- as.data.frame(res)
    df$year <- years
    df <- tidyr::pivot_longer(df, cols = -year, names_to = "area", values_to = "residual")
    fit_mat <- if (is.matrix(fit_val)) fit_val else matrix(fit_val, ncol = 1)
    df2 <- as.data.frame(fit_mat)
    df2$year <- years
    df2 <- tidyr::pivot_longer(df2, cols = -year, names_to = "area", values_to = "fitted_val")
    df <- merge(df, df2, by = c("year", "area"))
    facet_fn <- function(p) p + ggplot2::facet_wrap(~area)
  } else {
    df <- data.frame(year = years, residual = res, fitted_val = fit_val)
    facet_fn <- identity
  }

  # Remove NAs for plotting
  df <- df[is.finite(df$residual), ]

  # QQ plot
  if ("qq" %in% which || "all" %in% which) {
    p_qq <- ggplot2::ggplot(df, ggplot2::aes(sample = residual)) +
      ggplot2::stat_qq(colour = "#0072B2", size = 1.5) +
      ggplot2::stat_qq_line(colour = "grey40", linetype = "dashed") +
      ggplot2::labs(title = "Normal QQ Plot of Residuals",
                    x = "Theoretical Quantiles", y = "Sample Quantiles") +
      ggplot2::theme_bw()
    p_qq <- facet_fn(p_qq)
    plots$qq <- p_qq
    if (show) print(p_qq)
  }

  # Residuals vs fitted
  if ("fitted" %in% which || "all" %in% which) {
    p_rvf <- ggplot2::ggplot(df, ggplot2::aes(x = fitted_val, y = residual)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey") +
      ggplot2::geom_point(colour = "#D55E00", alpha = 0.7, size = 1.5) +
      ggplot2::geom_smooth(method = "loess", formula = y ~ x,
                            se = FALSE, colour = "#009E73", linewidth = 0.8) +
      ggplot2::labs(title = "Residuals vs Fitted CPUE",
                    x = "Fitted CPUE", y = "Log-Residual") +
      ggplot2::theme_bw()
    p_rvf <- facet_fn(p_rvf)
    plots$fitted <- p_rvf
    if (show) print(p_rvf)
  }

  # Histogram
  if ("histogram" %in% which || "all" %in% which) {
    p_hist <- ggplot2::ggplot(df, ggplot2::aes(x = residual)) +
      ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                               bins = max(10, nrow(df) %/% 3),
                               fill = "#56B4E9", colour = "white", alpha = 0.7) +
      ggplot2::stat_function(fun = dnorm,
                              args = list(mean = mean(df$residual),
                                          sd = sd(df$residual)),
                              colour = "grey30", linewidth = 0.8) +
      ggplot2::labs(title = "Histogram of Log-Residuals",
                    x = "Log-Residual", y = "Density") +
      ggplot2::theme_bw()
    p_hist <- facet_fn(p_hist)
    plots$histogram <- p_hist
    if (show) print(p_hist)
  }

  invisible(plots)
}


#' Plot Biomass Trajectory with Uncertainty
#'
#' Plots the estimated biomass trajectory with optional confidence bands
#' from \code{sdreport} standard errors, plus reference point lines.
#'
#' @param model_fit A fitted ProductionModel object.
#' @param ci Numeric, confidence level for bands (default: 0.95).
#'   Set to NULL to suppress bands.
#' @param show_ref Logical, overlay BMSY reference line (default: TRUE).
#' @param show Logical, whether to print the plot (default: TRUE).
#' @return A ggplot object (invisibly).
#' @export
plot_biomass <- function(model_fit, ci = 0.95, show_ref = TRUE, show = TRUE) {
  if (!inherits(model_fit, "ProductionModel")) stop("Input must be a ProductionModel object")
  if (!model_fit$fitted) stop("Model must be fitted before plotting")

  results <- model_fit$results
  years   <- model_fit$data$years
  bio     <- results$biomass
  bio_se  <- results$biomass_se

  if (is.matrix(bio)) {
    areas <- colnames(bio)
    df <- as.data.frame(bio)
    df$year <- years
    df <- tidyr::pivot_longer(df, cols = -year, names_to = "area", values_to = "biomass_val")
    if (!is.null(bio_se) && is.matrix(bio_se)) {
      se_df <- as.data.frame(bio_se)
      se_df$year <- years
      se_df <- tidyr::pivot_longer(se_df, cols = -year, names_to = "area", values_to = "se")
      df <- merge(df, se_df, by = c("year", "area"))
    }
  } else {
    df <- data.frame(year = years, biomass_val = bio)
    if (!is.null(bio_se)) df$se <- bio_se
  }

  # CI bands (log-normal approximation)
  if (!is.null(ci) && "se" %in% names(df)) {
    z <- qnorm(1 - (1 - ci) / 2)
    # Log-normal CI: exp(log(B) +/- z * se/B)
    df$lower <- df$biomass_val * exp(-z * df$se / df$biomass_val)
    df$upper <- df$biomass_val * exp( z * df$se / df$biomass_val)
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = year, y = biomass_val))

  if ("lower" %in% names(df)) {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                                   fill = "#0072B2", alpha = 0.2)
  }

  p <- p +
    ggplot2::geom_line(colour = "#0072B2", linewidth = 1.0) +
    ggplot2::geom_point(colour = "#0072B2", size = 1.5) +
    ggplot2::labs(title = "Estimated Biomass Trajectory",
                  y = "Biomass (tonnes)", x = "Year") +
    ggplot2::theme_bw()

  if (show_ref && !is.null(results$bmsy)) {
    p <- p + ggplot2::geom_hline(yintercept = results$bmsy,
                                  linetype = "dashed", colour = "#E69F00") +
      ggplot2::annotate("text", x = min(years), y = results$bmsy,
                         label = "BMSY", vjust = -0.5, hjust = 0,
                         colour = "#E69F00", size = 3.5)
  }

  if (is.matrix(bio)) {
    p <- p + ggplot2::facet_wrap(~area, scales = "free_y")
  }

  if (show) print(p)
  invisible(p)
}


# Declare global variables to appease R CMD check for NSE in ggplot2/tidyr
utils::globalVariables(c(
  "year", "area", "label", "observed", "fitted", "residuals", "harvest_rate",
  "residual", "fitted_val", "biomass_val", "se", "lower", "upper", "density"
))
