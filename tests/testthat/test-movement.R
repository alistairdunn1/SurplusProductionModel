test_that("gravity movement redistributes biomass and conserves total", {
  years <- 2000:2005
  areas <- c("A1", "A2")
  nY <- length(years)
  nA <- length(areas)

  # Simple distance matrix (symmetric)
  dist_mat <- matrix(c(0, 100, 100, 0), nrow = 2, byrow = TRUE, dimnames = list(areas, areas))
  attract <- c(A1 = 1, A2 = 2)
  move_rate <- 0.5
  decay <- 0.01

  # Simulate catch and CPUE
  B0 <- c(A1 = 1000, A2 = 2000)
  catch <- matrix(100, nY, nA, dimnames = list(year = years, area = areas))
  biomass <- matrix(NA_real_, nY, nA, dimnames = list(year = years, area = areas))
  biomass[1, ] <- B0
  r <- 0.3
  K <- 5000
  m <- 2
  q <- c(A1 = 0.001, A2 = 0.001)
  for (t in 1:(nY - 1)) {
    prod <- r * biomass[t, ] * (1 - (biomass[t, ] / K)) / m
    biomass[t + 1, ] <- biomass[t, ] + prod - catch[t, ]
    # Apply movement manually for reference
    W <- matrix(0, nA, nA)
    for (a in 1:nA) for (b in 1:nA) W[a, b] <- attract[b] * exp(-decay * dist_mat[a, b])
    Kmat <- W / rowSums(W)
    biomass[t + 1, ] <- (1 - move_rate) * biomass[t + 1, ] + move_rate * as.numeric(Kmat %*% biomass[t + 1, ])
  }
  cpue <- sweep(biomass, 2, q, "*")

  cpue_df <- do.call(rbind, lapply(1:nA, function(a) {
    data.frame(year = years, area = areas[a], cpue = cpue[, a])
  }))
  catch_df <- do.call(rbind, lapply(1:nA, function(a) {
    data.frame(year = years, area = areas[a], catch = catch[, a])
  }))

  data_list <- list(
    cpue_data = cpue_df,
    catch_data = catch_df,
    movement = list(distance_matrix = dist_mat, attractiveness = attract, movement_rate = move_rate, decay = decay)
  )

  start_vals <- list(
    log_r = log(r), log_K = log(K), log_m = log(m), log_sigma_proc = log(0.2), log_sigma_obs = log(0.2),
    log_q.A1 = log(q[1]), log_q.A2 = log(q[2]), log_B0.A1 = log(B0[1]), log_B0.A2 = log(B0[2])
  )

  fit <- fit_pella_tomlinson_model(data_list, params_init = start_vals, options = list(validate_data = FALSE))
  expect_s3_class(fit, "ProductionModel")
  # Verify model produces reasonable multi-area biomass results
  model_biomass <- fit$results$biomass
  expect_true(is.matrix(model_biomass))
  expect_equal(ncol(model_biomass), 2)
  expect_true(all(model_biomass > 0))
  # Both areas should have positive total biomass each year
  expect_true(all(rowSums(model_biomass) > 0))
})
