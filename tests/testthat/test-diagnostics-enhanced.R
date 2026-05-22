# Tests for enhanced diagnostic plots (plot_residuals, plot_biomass)
# and improved summary.ProductionModel

# ---------- Helper: fit a small model for testing -----------------------
fit_diag_model <- function() {
  set.seed(42)
  years <- 2010:2018
  true_r <- 0.3; true_K <- 5000; true_q <- 0.001; true_B0 <- 4000
  biomass <- numeric(length(years))
  biomass[1] <- true_B0
  catch_vals <- rep(600, length(years))
  for (t in seq_len(length(years) - 1)) {
    prod <- true_r * biomass[t] * (1 - biomass[t] / true_K)
    biomass[t + 1] <- max(100, biomass[t] + prod - catch_vals[t])
  }
  cpue <- true_q * biomass * exp(rnorm(length(years), 0, 0.05))

  data_list <- list(
    cpue_data  = data.frame(year = years, cpue = cpue),
    catch_data = data.frame(year = years, catch = catch_vals)
  )
  opts <- list(silent = TRUE, validate_data = FALSE,
               control = list(eval.max = 1000, iter.max = 500))
  fit_pella_tomlinson_model(data_list, options = opts)
}


# ========================  plot_residuals  ==============================

test_that("plot_residuals refuses unfitted or wrong class", {
  expect_error(plot_residuals("not a model"), "ProductionModel")
  unfitted <- structure(list(fitted = FALSE), class = "ProductionModel")
  expect_error(plot_residuals(unfitted), "fitted")
})

test_that("plot_residuals returns list of ggplots – all", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("ggplot2")
  model <- tryCatch(fit_diag_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  plots <- plot_residuals(model, which = "all", show = FALSE)
  expect_type(plots, "list")
  expect_true("qq" %in% names(plots))
  expect_true("fitted" %in% names(plots))
  expect_true("histogram" %in% names(plots))
  for (p in plots) expect_s3_class(p, "ggplot")
})

test_that("plot_residuals works for a single plot type", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("ggplot2")
  model <- tryCatch(fit_diag_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  qq_only <- plot_residuals(model, which = "qq", show = FALSE)
  expect_true("qq" %in% names(qq_only))
  expect_false("histogram" %in% names(qq_only))
})


# ========================  plot_biomass  ================================

test_that("plot_biomass refuses unfitted or wrong class", {
  expect_error(plot_biomass("not a model"), "ProductionModel")
  unfitted <- structure(list(fitted = FALSE), class = "ProductionModel")
  expect_error(plot_biomass(unfitted), "fitted")
})

test_that("plot_biomass returns a ggplot", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("ggplot2")
  model <- tryCatch(fit_diag_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  p <- plot_biomass(model, show = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_biomass without CI works", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("ggplot2")
  model <- tryCatch(fit_diag_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  p <- plot_biomass(model, ci = NULL, show = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_biomass with show_ref = FALSE works", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("ggplot2")
  model <- tryCatch(fit_diag_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  p <- plot_biomass(model, show_ref = FALSE, show = FALSE)
  expect_s3_class(p, "ggplot")
})


# ====================  summary.ProductionModel  =========================

test_that("summary.ProductionModel prints without error for fitted model", {
  skip_if_not_installed("RTMB")
  model <- tryCatch(fit_diag_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  out <- capture.output(res <- summary(model))
  expect_true(length(out) > 5)
  expect_true(any(grepl("Pella-Tomlinson", out)))
  expect_true(any(grepl("Reference Points", out)))
})

test_that("summary.ProductionModel can print depletion-based targets", {
  skip_if_not_installed("RTMB")
  model <- tryCatch(fit_diag_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  out <- capture.output(summary(model, biomass_target = 0.4, baseline = "K"))
  expect_true(any(grepl("User-Defined Biomass Targets", out)))
  expect_true(any(grepl("B_40%K", out, fixed = TRUE)))
  expect_true(any(grepl("F_40%K", out, fixed = TRUE)))
})

test_that("summary and print methods can use package-level defaults", {
  skip_if_not_installed("RTMB")
  model <- tryCatch(fit_diag_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  old_target <- getOption("SurplusProductionModel.biomass_target_default")
  old_baseline <- getOption("SurplusProductionModel.baseline_default")
  on.exit({
    do.call(options, setNames(list(old_target), "SurplusProductionModel.biomass_target_default"))
    do.call(options, setNames(list(old_baseline), "SurplusProductionModel.baseline_default"))
  }, add = TRUE)
  set_reference_point_defaults(biomass_target = 0.4, baseline = "K")

  out_summary <- capture.output(summary(model))
  out_print <- capture.output(print(model))

  expect_true(any(grepl("User-Defined Biomass Targets", out_summary)))
  expect_true(any(grepl("B_40%K", out_summary, fixed = TRUE)))
  expect_true(any(grepl("F_40%K", out_summary, fixed = TRUE)))

  expect_true(any(grepl("User-Defined Biomass Targets", out_print)))
  expect_true(any(grepl("B_40%K", out_print, fixed = TRUE)))
  expect_true(any(grepl("F_40%K", out_print, fixed = TRUE)))
})

test_that("summary returns invisible list for fitted model", {
  skip_if_not_installed("RTMB")
  model <- tryCatch(fit_diag_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  invisible_capture <- capture.output(res <- summary(model))
  expect_type(res, "list")
  expect_true("parameters" %in% names(res))
  expect_true("results" %in% names(res))
})

test_that("summary returns NULL for unfitted model", {
  unfitted <- structure(list(fitted = FALSE), class = "ProductionModel")
  expect_null(suppressMessages(summary(unfitted)))
})


# ====================  print.ProductionModel (enhanced) ==================

test_that("print.ProductionModel works for fitted model", {
  skip_if_not_installed("RTMB")
  model <- tryCatch(fit_diag_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  out <- capture.output(print(model))
  expect_true(any(grepl("Production Model Object", out)))
  expect_true(any(grepl("Convergence", out)))
  expect_true(any(grepl("Key Parameters", out)))
})

test_that("print.ProductionModel works for unfitted model", {
  unfitted <- structure(
    list(fitted = FALSE, model_type = "pella_tomlinson",
         creation_date = Sys.time(), data = list(), parameters = list(),
         results = list()),
    class = "ProductionModel"
  )
  out <- capture.output(print(unfitted))
  expect_true(any(grepl("Production Model Object", out)))
  expect_true(any(grepl("Fitted: No", out)))
})
