test_that("multi-index multi-area fit runs and produces plausible outputs", {
  skip_on_cran()
  set.seed(123)

  years <- 2000:2010
  areas <- c("A1", "A2")
  labels <- c("LL", "Trawl")
  nY <- length(years)
  nA <- length(areas)
  nL <- length(labels)

  # True parameters
  r <- 0.3
  K <- 5000
  m <- 2
  B0 <- c(A1 = 3000, A2 = 3500)
  q_true <- matrix(
    c(
      0.0004, 0.00025,
      0.0005, 0.0003
    ),
    nrow = 2, byrow = TRUE,
    dimnames = list(area = areas, label = labels)
  )

  # Simulate biomass per area (deterministic Schaefer dynamics)
  B <- matrix(NA_real_, nY, nA, dimnames = list(year = years, area = areas))
  B[1, ] <- B0
  C <- matrix(0, nY, nA, dimnames = list(year = years, area = areas))
  C[, "A1"] <- 200
  C[, "A2"] <- 240
  for (t in 1:(nY - 1)) {
    P <- r * B[t, ] * (1 - (B[t, ] / K)) # m=2
    B[t + 1, ] <- pmax(B[t, ] + P - C[t, ], 0.01)
  }

  # Simulate CPUE for each index with lognormal observation error
  sigma_obs <- 0.2
  cpue_arr <- array(NA_real_, dim = c(nY, nA, nL), dimnames = list(year = years, area = areas, label = labels))
  for (a in 1:nA) {
    for (l in 1:nL) {
      mu <- q_true[a, l] * B[, a]
      cpue_arr[, a, l] <- mu * exp(rnorm(nY, 0, sigma_obs))
    }
  }

  # Build input data.frames
  cpue_df <- do.call(rbind, lapply(1:nA, function(a) {
    do.call(rbind, lapply(1:nL, function(l) {
      data.frame(year = years, area = areas[a], label = labels[l], cpue = cpue_arr[, a, l])
    }))
  }))
  catch_df <- do.call(rbind, lapply(1:nA, function(a) {
    data.frame(year = years, area = areas[a], catch = C[, a])
  }))

  data_list <- list(cpue_data = cpue_df, catch_data = catch_df)

  # Starting values (let generator build per-index q and per-area B0)
  processed <- SurplusProductionModel:::preprocess_model_data(cpue_df, catch_df)
  start_vals <- SurplusProductionModel:::generate_starting_values(processed)

  # Ensure required names exist for per-index q
  expect_true(all(paste0("log_q.", rep(areas, each = nL), ".", rep(labels, times = nA)) %in% names(start_vals)))
  expect_true(all(paste0("log_B0.", areas) %in% names(start_vals)))

  # Fit
  fit <- fit_pella_tomlinson_model(data_list, params_init = start_vals, options = list(control = list(iter.max = 200, eval.max = 400), validate_data = FALSE))

  # Basic checks
  expect_s3_class(fit, "ProductionModel")
  expect_true(fit$fitted)
  expect_true(is.list(fit$results))
  expect_true(all(c("biomass", "fitted_cpue", "residuals", "harvest_rate") %in% names(fit$results)))

  # Shapes expected in results for multi-area
  expect_true(is.matrix(fit$results$biomass))
  expect_true(
    is.matrix(fit$results$fitted_cpue) ||
      (is.array(fit$results$fitted_cpue) && length(dim(fit$results$fitted_cpue)) == 3)
  )

  # Plausibility: biomass positive and not exploding
  expect_true(all(is.finite(fit$results$biomass)))
  expect_true(all(fit$results$biomass > 0))
})
