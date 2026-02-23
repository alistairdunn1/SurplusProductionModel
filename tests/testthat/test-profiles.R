# Tests for profile likelihood confidence intervals

# ---------- Helper: fit a small model for testing -----------------------
fit_profile_model <- function() {
  set.seed(42)
  years <- 2010:2020
  true_r <- 0.3; true_K <- 5000; true_q <- 0.001; true_B0 <- 4000
  biomass <- numeric(length(years))
  biomass[1] <- true_B0
  catch_vals <- rep(500, length(years))
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


# ====================  Input validation  ================================

test_that("profile_likelihood refuses unfitted or wrong class", {
  expect_error(profile_likelihood("not a model"), "ProductionModel")
  unfitted <- structure(list(fitted = FALSE), class = "ProductionModel")
  expect_error(profile_likelihood(unfitted), "fitted")
})

test_that("profile_likelihood errors on unknown parameter", {
  skip_if_not_installed("RTMB")
  model <- tryCatch(fit_profile_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })
  expect_error(profile_likelihood(model, parameters = "nonexistent"),
               "not found")
})


# ====================  Basic profiling  =================================

test_that("profile_likelihood returns correct structure for single param", {
  skip_if_not_installed("RTMB")
  model <- tryCatch(fit_profile_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  prof <- profile_likelihood(model, parameters = "r", n_points = 5)

  expect_s3_class(prof, "profile_likelihood")
  expect_true("profiles" %in% names(prof))
  expect_true("ci" %in% names(prof))
  expect_true("mle" %in% names(prof))
  expect_true("r" %in% names(prof$profiles))

  # Profile data.frame structure
  df <- prof$profiles$r
  expect_true(is.data.frame(df))
  expect_true(all(c("value", "nll", "delta_nll", "converged") %in% names(df)))
  expect_true(nrow(df) >= 10)  # at least 2*5+1 points

  # MLE should be positive

  expect_true(prof$mle[["r"]] > 0)

  # CI should be a length-2 vector
  expect_length(prof$ci$r, 2)
  expect_true(all(names(prof$ci$r) %in% c("lower", "upper")))
})

test_that("profile_likelihood MLE has delta_nll ~ 0", {
  skip_if_not_installed("RTMB")
  model <- tryCatch(fit_profile_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  prof <- profile_likelihood(model, parameters = "r", n_points = 5)
  df <- prof$profiles$r

  # At the MLE point, delta_nll should be very close to 0
  mle_row <- which.min(df$delta_nll)
  expect_true(df$delta_nll[mle_row] < 0.01)
})

test_that("profile CI brackets the MLE", {
  skip_if_not_installed("RTMB")
  model <- tryCatch(fit_profile_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  prof <- profile_likelihood(model, parameters = "r", n_points = 10)
  ci <- prof$ci$r
  mle <- prof$mle[["r"]]

  # If both bounds were found, MLE should be within the CI
  if (!is.na(ci["lower"])) expect_true(ci["lower"] < mle)
  if (!is.na(ci["upper"])) expect_true(ci["upper"] > mle)
})


# ====================  Multiple parameters  =============================

test_that("profile_likelihood works with multiple parameters", {
  skip_if_not_installed("RTMB")
  model <- tryCatch(fit_profile_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  prof <- profile_likelihood(model, parameters = c("r", "K"), n_points = 5)

  expect_true("r" %in% names(prof$profiles))
  expect_true("K" %in% names(prof$profiles))
  expect_length(prof$ci, 2)
  expect_length(prof$mle, 2)
})


# ====================  Derived quantities  ==============================

test_that("profile_likelihood works for derived quantity MSY", {
  skip_if_not_installed("RTMB")
  model <- tryCatch(fit_profile_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  prof <- profile_likelihood(model, parameters = "MSY", n_points = 5)

  expect_true("MSY" %in% names(prof$profiles))
  df <- prof$profiles$MSY
  expect_true(is.data.frame(df))
  expect_true(nrow(df) >= 5)
  expect_true(prof$mle[["MSY"]] > 0)
})

test_that("profile_likelihood works for BMSY", {
  skip_if_not_installed("RTMB")
  model <- tryCatch(fit_profile_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  prof <- profile_likelihood(model, parameters = "BMSY", n_points = 5)
  expect_true("BMSY" %in% names(prof$profiles))
  expect_true(prof$mle[["BMSY"]] > 0)
})


# ====================  Print and plot methods  ==========================

test_that("print.profile_likelihood doesn't error", {
  skip_if_not_installed("RTMB")
  model <- tryCatch(fit_profile_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  prof <- profile_likelihood(model, parameters = "r", n_points = 5)
  expect_output(print(prof), "Profile Likelihood")
  expect_output(print(prof), "MLE")
})

test_that("plot.profile_likelihood returns a ggplot", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("ggplot2")
  model <- tryCatch(fit_profile_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  prof <- profile_likelihood(model, parameters = c("r", "K"), n_points = 5)
  p <- plot(prof)
  expect_s3_class(p, "ggplot")
})


# ====================  CI level parameter  ==============================

test_that("profile CI widens with higher confidence level", {
  skip_if_not_installed("RTMB")
  model <- tryCatch(fit_profile_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  prof90 <- profile_likelihood(model, parameters = "r", ci_level = 0.90, n_points = 10)
  prof99 <- profile_likelihood(model, parameters = "r", ci_level = 0.99, n_points = 10)

  ci90 <- prof90$ci$r
  ci99 <- prof99$ci$r

  # 99% CI should be at least as wide as 90% CI (where both bounds exist)
  if (!is.na(ci90["lower"]) && !is.na(ci99["lower"])) {
    expect_true(ci99["lower"] <= ci90["lower"])
  }
  if (!is.na(ci90["upper"]) && !is.na(ci99["upper"])) {
    expect_true(ci99["upper"] >= ci90["upper"])
  }
})
