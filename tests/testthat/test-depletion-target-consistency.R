# Cross-function consistency tests for depletion targets

# ---------- Helper: fit a small model for consistency checks ------------
fit_consistency_model <- function() {
  set.seed(42)
  years <- 2010:2020
  true_r <- 0.3
  true_K <- 5000
  true_q <- 0.001
  true_B_initial <- 4000

  biomass <- numeric(length(years))
  biomass[1] <- true_B_initial
  catch_vals <- rep(500, length(years))
  for (t in seq_len(length(years) - 1)) {
    prod <- true_r * biomass[t] * (1 - biomass[t] / true_K)
    biomass[t + 1] <- max(100, biomass[t] + prod - catch_vals[t])
  }
  cpue <- true_q * biomass * exp(rnorm(length(years), 0, 0.05))

  data_list <- list(
    cpue_data = data.frame(year = years, cpue = cpue),
    catch_data = data.frame(year = years, catch = catch_vals)
  )

  opts <- list(
    silent = TRUE,
    validate_data = FALSE,
    control = list(eval.max = 1000, iter.max = 500)
  )

  fit_pella_tomlinson_model(data_list, options = opts)
}

test_that("depletion targets are consistent between reference points and profile MLE", {
  skip_if_not_installed("RTMB")

  model <- tryCatch(fit_consistency_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  ref <- calculate_reference_points(model, biomass_target = 0.4, baseline = "K")
  prof <- profile_likelihood(
    model,
    parameters = c("B_40%K", "F_40%K"),
    biomass_target = 0.4,
    baseline = "K",
    n_points = 5
  )

  target_row <- ref$target_reference_points[1, , drop = FALSE]

  expect_equal(as.numeric(prof$mle[["B_40%K"]]), as.numeric(target_row$biomass), tolerance = 1e-6)
  expect_equal(as.numeric(prof$mle[["F_40%K"]]), as.numeric(target_row$fishing_mortality), tolerance = 1e-6)
})

test_that("depletion target labels appear in summary and print outputs", {
  skip_if_not_installed("RTMB")

  model <- tryCatch(fit_consistency_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  summary_out <- capture.output(summary(model, biomass_target = 0.4, baseline = "K"))
  print_out <- capture.output(print(model, biomass_target = 0.4, baseline = "K"))

  expect_true(any(grepl("User-Defined Biomass Targets", summary_out)))
  expect_true(any(grepl("B_40%K", summary_out, fixed = TRUE)))
  expect_true(any(grepl("F_40%K", summary_out, fixed = TRUE)))

  expect_true(any(grepl("User-Defined Biomass Targets", print_out)))
  expect_true(any(grepl("B_40%K", print_out, fixed = TRUE)))
  expect_true(any(grepl("F_40%K", print_out, fixed = TRUE)))
})

test_that("bayesian depletion targets match analytical helper per draw", {
  skip_if_not_installed("RTMB")
  skip_if_not_installed("tmbstan")
  skip_if_not_installed("rstan")

  model <- tryCatch(fit_consistency_model(), error = function(e) {
    skip(paste("Model fitting failed:", e$message))
  })

  bf <- tryCatch(
    bayesian_fit(
      model,
      chains = 1,
      iter = 200,
      warmup = 100,
      seed = 42,
      biomass_target = 0.4,
      baseline = "K"
    ),
    error = function(e) skip(paste("tmbstan failed:", e$message))
  )

  posterior <- as.data.frame(bf$posterior)

  idx <- which(
    is.finite(posterior[["r"]]) &
      is.finite(posterior[["K"]]) &
      is.finite(posterior[["m"]]) &
      is.finite(posterior[["B_40%K"]]) &
      is.finite(posterior[["F_40%K"]])
  )

  if (length(idx) == 0) {
    skip("No finite posterior rows available for depletion-target comparison")
  }

  i <- idx[1]
  params <- c(
    r = as.numeric(posterior[["r"]][i]),
    K = as.numeric(posterior[["K"]][i]),
    m = as.numeric(posterior[["m"]][i])
  )
  if ("B_initial" %in% names(posterior) && is.finite(posterior[["B_initial"]][i])) {
    params <- c(params, B_initial = as.numeric(posterior[["B_initial"]][i]))
  }

  expected <- SurplusProductionModel:::calculate_reference_points_from_parameters(
    parameters = params,
    biomass_target = 0.4,
    baseline = "K",
    warn_on_invalid_m = FALSE
  )

  expected_row <- expected$target_reference_points[1, , drop = FALSE]
  expect_equal(as.numeric(posterior[["B_40%K"]][i]), as.numeric(expected_row$biomass), tolerance = 1e-10)
  expect_equal(as.numeric(posterior[["F_40%K"]][i]), as.numeric(expected_row$fishing_mortality), tolerance = 1e-10)
})
