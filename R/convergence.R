#' Convergence Diagnostics for Surplus Production Model
#'
#' Functions for assessing model convergence reliability, including
#' jitter tests and retrospective analysis.
#'
#' @name convergence
NULL

#' Jitter Test for Convergence Reliability
#'
#' Re-fits the model from multiple randomly perturbed starting values to
#' assess whether the optimizer consistently finds the same minimum.
#' Large variation in final negative log-likelihood or parameter estimates
#' across jittered starts indicates convergence problems.
#'
#' @param model_fit A fitted \code{ProductionModel} object.
#' @param n_jitter Integer, number of jittered re-fits (default: 10).
#' @param jitter_sd Numeric, standard deviation of the normal perturbation
#'   added to log-scale starting values (default: 0.2).
#' @param tolerance Numeric, maximum acceptable difference in negative
#'   log-likelihood from the best run before flagging a failure
#'   (default: 2.0, roughly one likelihood unit).
#'
#' @return A list with class \code{"jitter_test"} containing:
#'   \describe{
#'     \item{nll}{Numeric vector of negative log-likelihoods for each run.}
#'     \item{convergence}{Integer vector of convergence codes (0 = success).}
#'     \item{parameters}{Data frame of fitted parameter values (rows = runs,
#'       columns = parameters).}
#'     \item{best_nll}{Numeric, the smallest NLL achieved.}
#'     \item{n_within_tol}{Integer, number of runs within \code{tolerance}
#'       of the best NLL.}
#'     \item{proportion_converged}{Numeric, fraction of runs that converged
#'       and fell within tolerance.}
#'     \item{passed}{Logical, TRUE if > 50\% of runs converge to the same
#'       minimum (within tolerance).}
#'   }
#'
#' @details
#' Each jitter run creates starting values by adding
#' \eqn{N(0, \sigma_j)} to every log-scale parameter, then re-runs the
#' full optimisation pipeline.  The original fitted values are included
#' as run 1 (un-jittered).
#'
#' A model that consistently returns the same NLL (within tolerance)
#' from diverse starting points is more trustworthy than one that
#' finds different minima.
#'
#' @examples
#' \dontrun{
#' jt <- jitter_test(fitted_model, n_jitter = 20)
#' print(jt)
#' plot(jt)
#' }
#'
#' @export
jitter_test <- function(model_fit, n_jitter = 10, jitter_sd = 0.2,
                        tolerance = 2.0) {
  if (!inherits(model_fit, "ProductionModel")) {
    stop("Input must be a fitted ProductionModel object")
  }
  if (!model_fit$fitted) stop("Model must be fitted before jitter testing")

  obj <- model_fit$results$rtmb_obj
  if (is.null(obj)) stop("rtmb_obj not found in results -- was the model fitted with RTMB?")

  base_par  <- model_fit$results$opt_par
  if (is.null(base_par)) {
    stop("Optimised parameters not found in results")
  }
  ctrl      <- list(eval.max = 1000, iter.max = 500)
  n_total   <- n_jitter + 1L   # include the original

  nll_vec  <- numeric(n_total)
  conv_vec <- integer(n_total)
  par_mat  <- matrix(NA_real_, nrow = n_total, ncol = length(base_par))
  colnames(par_mat) <- names(base_par)

  for (i in seq_len(n_total)) {
    if (i == 1L) {
      # Run 1 represents the original fit; do not optimise it again because
      # doing so can move to a different minimum and no longer reports the
      # likelihood or convergence status of the supplied model.
      nll_vec[i] <- model_fit$results$likelihood
      conv_vec[i] <- model_fit$results$convergence
      par_mat[i, ] <- base_par
      next
    }

    start_par <- base_par + rnorm(length(base_par), 0, jitter_sd)
    opt <- tryCatch(
      nlminb(start_par, obj$fn, obj$gr, control = ctrl),
      error = function(e) list(objective = NA_real_, convergence = 99L, par = rep(NA_real_, length(base_par)))
    )
    nll_vec[i]   <- opt$objective
    conv_vec[i]  <- opt$convergence
    par_mat[i, ] <- opt$par
  }

  best_nll      <- min(nll_vec, na.rm = TRUE)
  within_tol    <- which(!is.na(nll_vec) & conv_vec == 0 &
                           abs(nll_vec - best_nll) <= tolerance)
  n_within      <- length(within_tol)
  prop_converged <- n_within / n_total

  out <- list(
    nll                  = nll_vec,
    convergence          = conv_vec,
    parameters           = as.data.frame(par_mat),
    best_nll             = best_nll,
    n_within_tol         = n_within,
    n_total              = n_total,
    proportion_converged = prop_converged,
    tolerance            = tolerance,
    jitter_sd            = jitter_sd,
    passed               = prop_converged > 0.5
  )
  class(out) <- "jitter_test"
  out
}


#' @export
print.jitter_test <- function(x, ...) {
  cat("Jitter Test Results\n")
  cat("===================\n")
  cat("Runs:", x$n_total, "(1 original +", x$n_total - 1L, "jittered)\n")
  cat("Jitter SD:", x$jitter_sd, "\n")
  cat("Best NLL:", sprintf("%.4f", x$best_nll), "\n")
  cat("Runs within tolerance (", x$tolerance, "):",
      x$n_within_tol, "/", x$n_total,
      sprintf(" (%.0f%%)\n", x$proportion_converged * 100))
  cat("NLL range:", sprintf("%.4f - %.4f", min(x$nll, na.rm = TRUE),
                             max(x$nll, na.rm = TRUE)), "\n")
  cat("Result:", if (x$passed) "PASSED" else "FAILED", "\n")
  invisible(x)
}


#' Plot Jitter Test Results
#'
#' @param x A \code{jitter_test} object.
#' @param ... Additional arguments (not used).
#' @return A ggplot object (invisibly).
#' @export
plot.jitter_test <- function(x, ...) {
  df <- data.frame(
    run = seq_along(x$nll),
    nll = x$nll,
    converged = x$convergence == 0
  )
  p <- ggplot2::ggplot(df, ggplot2::aes(x = run, y = nll)) +
    ggplot2::geom_point(ggplot2::aes(colour = converged), size = 2.5) +
    ggplot2::geom_hline(yintercept = x$best_nll, linetype = "dashed",
                         colour = "grey40") +
    ggplot2::geom_hline(yintercept = x$best_nll + x$tolerance,
                         linetype = "dotted", colour = "red") +
    ggplot2::scale_colour_manual(values = c("TRUE" = "#009E73", "FALSE" = "#D55E00"),
                                  labels = c("TRUE" = "Converged", "FALSE" = "Failed")) +
    ggplot2::labs(title = "Jitter Test: NLL by Start",
                  x = "Run", y = "Negative Log-Likelihood",
                  colour = "Status") +
    ggplot2::theme_bw()
  print(p)
  invisible(p)
}


#' Retrospective Analysis
#'
#' Performs a retrospective analysis by sequentially removing the most recent
#' year(s) of data and re-fitting the model.  This reveals systematic bias
#' in parameter estimates or biomass trajectories caused by incoming data
#' (Mohn's rho).
#'
#' @param model_fit A fitted \code{ProductionModel} object.
#' @param n_peels Integer, number of years to peel (default: 5).
#' @param options List of options passed to \code{fit_pella_tomlinson_model()}
#'   for each re-fit (e.g. \code{fixed_params}, \code{control}).  If NULL,
#'   defaults from the original fit are inferred.
#'
#' @return A list with class \code{"retro_analysis"} containing:
#'   \describe{
#'     \item{peels}{List of fitted \code{ProductionModel} objects, one per peel
#'       (element 1 = full data, element 2 = minus 1 year, etc.).}
#'     \item{terminal_biomass}{Named numeric vector of terminal-year biomass
#'       estimates for each peel.}
#'     \item{mohns_rho}{Mohn's rho for terminal biomass. Values far from zero
#'       indicate retrospective bias.}
#'   }
#'
#' @details
#' Mohn's rho is calculated as the mean relative difference between the
#' peeled terminal biomass and the corresponding biomass from the full model:
#' \deqn{\rho = \frac{1}{n_{\text{peels}}} \sum_{p=1}^{n_{\text{peels}}}
#' \frac{B^{(-p)}_{T-p} - B^{(0)}_{T-p}}{B^{(0)}_{T-p}}}
#'
#' Values of \eqn{|\rho| > 0.20} are generally considered problematic.
#'
#' @examples
#' \dontrun{
#' retro <- retrospective_analysis(fitted_model, n_peels = 5)
#' print(retro)
#' plot(retro)
#' }
#'
#' @export
retrospective_analysis <- function(model_fit, n_peels = 5, options = NULL) {
  if (!inherits(model_fit, "ProductionModel")) {
    stop("Input must be a fitted ProductionModel object")
  }
  if (!model_fit$fitted) stop("Model must be fitted before retrospective analysis")

  # Extract original data from the fitted model
  orig_data <- model_fit$data
  all_years <- orig_data$years
  n_years   <- length(all_years)

  if (n_peels >= n_years - 3) {
    stop("n_peels (", n_peels, ") too large for a time series of ", n_years, " years")
  }

  # Reconstruct the data list that was passed to fit_pella_tomlinson_model()
  # We need cpue_data and catch_data as data.frames
  areas  <- orig_data$areas
  cpue   <- orig_data$cpue
  ctch   <- orig_data$catch
  has_labels <- is.array(cpue) && length(dim(cpue)) == 3

  # Build cpue_data data.frame
  if (has_labels) {
    cpue_df <- do.call(rbind, lapply(seq_along(areas), function(ia) {
      do.call(rbind, lapply(dimnames(cpue)[[3]], function(l) {
        data.frame(year = all_years, area = areas[ia],
                   label = l, cpue = cpue[, ia, l],
                   stringsAsFactors = FALSE)
      }))
    }))
  } else if (is.matrix(cpue)) {
    cpue_df <- do.call(rbind, lapply(seq_along(areas), function(ia) {
      data.frame(year = all_years, area = areas[ia],
                 cpue = cpue[, ia], stringsAsFactors = FALSE)
    }))
  } else {
    cpue_df <- data.frame(year = all_years, cpue = cpue, stringsAsFactors = FALSE)
  }

  # Build catch_data data.frame
  if (is.matrix(ctch)) {
    catch_df <- do.call(rbind, lapply(seq_along(areas), function(ia) {
      data.frame(year = all_years, area = areas[ia],
                 catch = ctch[, ia], stringsAsFactors = FALSE)
    }))
  } else {
    catch_df <- data.frame(year = all_years, catch = ctch, stringsAsFactors = FALSE)
  }

  # Default options for re-fits: silent, no validation (already validated)
  if (is.null(options)) {
    options <- list(silent = TRUE, validate_data = FALSE,
                    control = list(eval.max = 1000, iter.max = 500))
  }

  # Peel 0 = full model (already fitted)
  peels <- vector("list", n_peels + 1)
  peels[[1]] <- model_fit

  for (p in seq_len(n_peels)) {
    max_year <- max(all_years) - p
    cpue_p  <- cpue_df[cpue_df$year <= max_year, , drop = FALSE]
    catch_p <- catch_df[catch_df$year <= max_year, , drop = FALSE]

    data_p <- list(cpue_data = cpue_p, catch_data = catch_p)
    # Carry over movement if present
    if (!is.null(orig_data$movement_rate) || !is.null(model_fit$results$process_noise)) {
      # movement info not stored directly in data -- rebuild from original
    }

    peels[[p + 1]] <- tryCatch(
      fit_pella_tomlinson_model(data_p, options = options),
      error = function(e) {
        warning("Retrospective peel -", p, " failed: ", e$message)
        NULL
      }
    )
  }

  # Extract terminal biomass for each peel
  full_biomass <- model_fit$results$biomass
  terminal_bio <- numeric(n_peels + 1)
  names(terminal_bio) <- paste0("-", 0:n_peels)

  # Full model terminal biomass
  if (is.matrix(full_biomass)) {
    terminal_bio[1] <- sum(full_biomass[n_years, ])
  } else {
    terminal_bio[1] <- full_biomass[n_years]
  }

  # Peeled terminal biomasses and Mohn's rho calculation
  rho_terms <- numeric(n_peels)
  for (p in seq_len(n_peels)) {
    pfit <- peels[[p + 1]]
    if (is.null(pfit)) {
      terminal_bio[p + 1] <- NA
      rho_terms[p] <- NA
      next
    }
    peel_bio <- pfit$results$biomass
    peel_n   <- length(pfit$data$years)
    if (is.matrix(peel_bio)) {
      terminal_bio[p + 1] <- sum(peel_bio[peel_n, ])
    } else {
      terminal_bio[p + 1] <- peel_bio[peel_n]
    }

    # Full model biomass at the same year
    year_idx <- which(all_years == max(all_years) - p)
    full_at_year <- if (is.matrix(full_biomass)) sum(full_biomass[year_idx, ]) else full_biomass[year_idx]
    rho_terms[p] <- (terminal_bio[p + 1] - full_at_year) / full_at_year
  }

  mohns_rho <- mean(rho_terms, na.rm = TRUE)

  out <- list(
    peels             = peels,
    terminal_biomass  = terminal_bio,
    mohns_rho         = mohns_rho,
    n_peels           = n_peels,
    years             = all_years
  )
  class(out) <- "retro_analysis"
  out
}


#' @export
print.retro_analysis <- function(x, ...) {
  cat("Retrospective Analysis\n")
  cat("======================\n")
  cat("Peels:", x$n_peels, "\n")
  cat("Years:", min(x$years), "-", max(x$years), "\n\n")
  cat("Terminal biomass by peel:\n")
  for (i in seq_along(x$terminal_biomass)) {
    label <- names(x$terminal_biomass)[i]
    val   <- x$terminal_biomass[i]
    cat(sprintf("  Peel %s: %.1f\n", label, val))
  }
  cat(sprintf("\nMohn's rho: %.4f", x$mohns_rho))
  if (abs(x$mohns_rho) > 0.20) {
    cat(" ** WARNING: |rho| > 0.20 indicates retrospective bias **")
  }
  cat("\n")
  invisible(x)
}


#' Plot Retrospective Analysis
#'
#' Plots overlaid biomass trajectories from each retrospective peel.
#'
#' @param x A \code{retro_analysis} object.
#' @param ... Additional arguments (not used).
#' @return A ggplot object (invisibly).
#' @export
plot.retro_analysis <- function(x, ...) {
  # Build a long data.frame of biomass from each peel
  dfs <- list()
  for (i in seq_along(x$peels)) {
    pfit <- x$peels[[i]]
    if (is.null(pfit)) next
    bio <- pfit$results$biomass
    yrs <- pfit$data$years
    if (is.matrix(bio)) {
      total_bio <- rowSums(bio)
    } else {
      total_bio <- bio
    }
    dfs[[length(dfs) + 1]] <- data.frame(
      year    = yrs,
      biomass = total_bio,
      peel    = paste0("-", i - 1L),
      stringsAsFactors = FALSE
    )
  }
  df <- do.call(rbind, dfs)
  df$peel <- factor(df$peel, levels = paste0("-", 0:x$n_peels))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = year, y = biomass, colour = peel)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::labs(title = sprintf("Retrospective Analysis (Mohn's rho = %.3f)", x$mohns_rho),
                  x = "Year", y = "Total Biomass",
                  colour = "Peel") +
    ggplot2::theme_bw()
  print(p)
  invisible(p)
}

# Declare ggplot NSE variables for R CMD check
utils::globalVariables(c("run", "nll", "converged", "peel", "biomass"))
