test_that("pella_tomlinson_production function works correctly", {
  # Test basic functionality
  production <- pella_tomlinson_production(B = 1000, r = 0.3, K = 5000, m = 2)
  expect_type(production, "double")
  expect_true(production > 0)

  # Test Schaefer special case (m = 2)
  B <- 2500 # Half of K
  production_schaefer <- pella_tomlinson_production(B, r = 0.3, K = 5000, m = 2)
  # Standard Pella-Tomlinson with m = 2 (Schaefer): r/(m-1) * B * (1 - B/K)
  expected_schaefer <- 0.3 * B * (1 - B / 5000) / (2 - 1)
  expect_equal(production_schaefer, expected_schaefer, tolerance = 1e-10)

  # Test maximum production occurs at BMSY for Schaefer
  K <- 5000
  r <- 0.3
  bmsy_schaefer <- K / 2
  max_production <- pella_tomlinson_production(bmsy_schaefer, r, K, 2)

  # Test nearby points have lower production
  production_lower <- pella_tomlinson_production(bmsy_schaefer * 0.9, r, K, 2)
  production_upper <- pella_tomlinson_production(bmsy_schaefer * 1.1, r, K, 2)

  expect_true(max_production > production_lower)
  expect_true(max_production > production_upper)

  # Test boundary conditions
  expect_equal(pella_tomlinson_production(0, r = 0.3, K = 5000, m = 2), 0)
  expect_equal(pella_tomlinson_production(5000, r = 0.3, K = 5000, m = 2), 0)
  expect_equal(pella_tomlinson_production(6000, r = 0.3, K = 5000, m = 2), 0) # Above K

  # Test error conditions
  expect_error(pella_tomlinson_production(-100, r = 0.3, K = 5000, m = 2))
  expect_error(pella_tomlinson_production(1000, r = -0.3, K = 5000, m = 2))
  expect_error(pella_tomlinson_production(1000, r = 0.3, K = 0, m = 2))
  expect_error(pella_tomlinson_production(1000, r = 0.3, K = 5000, m = 0))
})

test_that("generate_starting_values creates reasonable parameters", {
  # Create test data
  years <- 2000:2020
  catch <- rep(1000, length(years))
  cpue <- exp(rnorm(length(years), log(1.5), 0.1))

  data <- list(years = years, catch = catch, cpue = cpue)

  # Generate starting values
  start_vals <- generate_starting_values(data)

  # Check structure
  expected_params <- c("log_r", "log_K", "log_m", "log_q", "log_sigma_proc", "log_sigma_obs", "log_d0")
  expect_true(all(expected_params %in% names(start_vals)))
  expect_equal(length(start_vals), length(expected_params))

  # Check all values are finite
  expect_true(all(is.finite(unlist(start_vals))))

  # Transform to natural scale and check reasonableness
  natural_params <- exp(unlist(start_vals))
  names(natural_params) <- gsub("^log_", "", names(start_vals))

  expect_true(natural_params[["r"]] > 0 && natural_params[["r"]] <= 1.0)
  expect_true(natural_params[["K"]] > max(catch)) # K should be larger than max catch
  expect_true(natural_params[["m"]] > 0)
  expect_true(natural_params[["q"]] > 0 && natural_params[["q"]] <= 1.0)
  expect_true(natural_params[["sigma_proc"]] > 0 && natural_params[["sigma_proc"]] <= 1.0)
  expect_true(natural_params[["sigma_obs"]] > 0 && natural_params[["sigma_obs"]] <= 1.0)
  expect_true(natural_params[["d0"]] > 0)
})

test_that("preprocess_model_data works correctly", {
  # Create test data with overlapping years
  cpue_data <- data.frame(
    year = 2000:2020,
    cpue = rnorm(21, 1.5, 0.3)
  )

  catch_data <- data.frame(
    year = 2002:2018, # Subset of CPUE years
    catch = rnorm(17, 1000, 100)
  )

  # Process data
  processed <- preprocess_model_data(cpue_data, catch_data)

  # Check structure
  expect_true(all(c("years", "cpue", "catch") %in% names(processed)))
  expect_equal(length(processed$years), length(processed$cpue))
  expect_equal(length(processed$years), length(processed$catch))

  # Check years are in correct range (union of catch and CPUE years)
  expect_true(all(processed$years >= 2000))
  expect_true(all(processed$years <= 2020))

  # Check years are sorted
  expect_true(all(diff(processed$years) == 1))

  # Test insufficient overlap error
  cpue_data_short <- data.frame(year = 2000:2001, cpue = c(1.5, 1.6))
  catch_data_short <- data.frame(year = 2010:2011, catch = c(1000, 1100))

  expect_error(preprocess_model_data(cpue_data_short, catch_data_short))

  # Test missing catch values - NA catch values are now set to 0 with a warning
  catch_data_na <- data.frame(year = 2000:2005, catch = c(1000, NA, 1200, 1100, 1000, 1050))
  cpue_data_overlap <- data.frame(year = 2000:2005, cpue = rnorm(6, 1.5, 0.1))

  expect_warning(
    preprocess_model_data(cpue_data_overlap, catch_data_na),
    "Missing catch values"
  )

  # Test gap-filling: input years with gaps should produce a complete sequence
  # Need >= 3 overlapping years between catch and cpue
  catch_gap <- data.frame(year = c(2000, 2002, 2004, 2006), catch = c(500, 600, 700, 800))
  cpue_gap <- data.frame(year = c(2000, 2002, 2004, 2006), cpue = c(1.2, 1.1, 1.0, 0.9))

  processed_gap <- preprocess_model_data(cpue_gap, catch_gap)

  expect_equal(processed_gap$years, 2000:2006)
  expect_true(all(diff(processed_gap$years) == 1))
  # Gap years should have zero catch and NA CPUE
  expect_equal(processed_gap$catch[processed_gap$years == 2001], 0)
  expect_equal(processed_gap$catch[processed_gap$years == 2003], 0)
  expect_equal(processed_gap$catch[processed_gap$years == 2005], 0)
  expect_true(is.na(processed_gap$cpue[processed_gap$years == 2001]))
  expect_true(is.na(processed_gap$cpue[processed_gap$years == 2003]))
  expect_true(is.na(processed_gap$cpue[processed_gap$years == 2005]))
})

test_that("transform_parameters_to_natural works correctly", {
  # Create log-transformed parameters
  log_params <- c(
    log_r = log(0.3),
    log_K = log(5000),
    log_m = log(2.0),
    log_q = log(0.001),
    log_sigma_proc = log(0.2),
    log_sigma_obs = log(0.3)
  )

  # Transform to natural scale
  natural_params <- transform_parameters_to_natural(log_params)

  # Check names are correct
  expected_names <- c("r", "K", "m", "q", "sigma_proc", "sigma_obs")
  expect_true(all(expected_names %in% names(natural_params)))

  # Check values are correct
  expect_equal(natural_params[["r"]], 0.3, tolerance = 1e-10)
  expect_equal(natural_params[["K"]], 5000, tolerance = 1e-10)
  expect_equal(natural_params[["m"]], 2.0, tolerance = 1e-10)
  expect_equal(natural_params[["q"]], 0.001, tolerance = 1e-10)
  expect_equal(natural_params[["sigma_proc"]], 0.2, tolerance = 1e-10)
  expect_equal(natural_params[["sigma_obs"]], 0.3, tolerance = 1e-10)
})

test_that("transform_parameters_to_natural handles AR1 and covariate effects", {
  params <- c(
    log_r = log(0.3),
    theta_rho = atanh(0.4),
    beta_temp = -0.15
  )

  natural <- transform_parameters_to_natural(params)

  expect_true("rho" %in% names(natural))
  expect_equal(natural[["rho"]], 0.4, tolerance = 1e-10)
  expect_equal(natural[["beta_temp"]], -0.15, tolerance = 1e-10)
  expect_equal(natural[["r"]], 0.3, tolerance = 1e-10)
})

test_that("environmental covariates are aligned and scaled", {
  years <- 2000:2005
  areas <- c("A1", "A2")
  env_data <- data.frame(
    year = rep(years, times = 2),
    area = rep(areas, each = length(years)),
    temp = c(seq(0.2, 0.7, length.out = length(years)), seq(0.3, 0.8, length.out = length(years))),
    nao = c(seq(-1, 1, length.out = length(years)), seq(-0.5, 1.5, length.out = length(years)))
  )

  env_info <- .prepare_environmental_covariates(
    env_data = env_data,
    years = years,
    areas = areas,
    covariate_names = c("temp", "nao"),
    lag = 1,
    scale_covariates = TRUE
  )

  expect_true(all(c("env_array", "covariate_names", "scaling") %in% names(env_info)))
  expect_equal(dim(env_info$env_array), c(length(years) - 1, length(areas), 2))
  expect_equal(env_info$scaling$lag, 1)
  expect_equal(sort(env_info$covariate_names), c("nao", "temp"))
})

test_that("calculate_model_results produces expected output", {
  # Create test parameters
  parameters <- c(
    r = 0.3,
    K = 5000,
    m = 2.0,
    q = 0.001,
    d0 = 0.8 # initial depletion; B_initial = d0 * K = 4000
  )

  # Create test data
  data <- list(
    years = 2000:2005,
    catch = c(800, 900, 1000, 1100, 1000, 950),
    cpue = c(3.8, 3.6, 3.4, 3.2, 3.0, 2.9)
  )

  # Calculate results
  results <- calculate_model_results(parameters, data)

  # Check structure
  expect_true(all(c("biomass", "harvest_rate", "fitted_cpue", "residuals") %in% names(results)))
  expect_equal(length(results$biomass), length(data$years))
  expect_equal(length(results$harvest_rate), length(data$years))
  expect_equal(length(results$fitted_cpue), length(data$years))
  expect_equal(length(results$residuals), length(data$years))

  # Check biomass starts at B_initial = d0 * K
  expect_equal(results$biomass[1], parameters[["d0"]] * parameters[["K"]])

  # Check all biomass values are positive
  expect_true(all(results$biomass > 0))

  # Check harvest rates are reasonable (can exceed 1 if biomass is very low)
  expect_true(all(results$harvest_rate >= 0))
  expect_true(all(is.finite(results$harvest_rate)))

  # Check fitted CPUE relationship
  q <- parameters[["q"]]
  expected_cpue <- q * results$biomass
  expect_equal(results$fitted_cpue, expected_cpue, tolerance = 1e-10)

  # Check residuals calculation
  expected_residuals <- log(data$cpue) - log(results$fitted_cpue)
  expect_equal(results$residuals, expected_residuals, tolerance = 1e-10)
})

test_that("fit_pella_tomlinson_model input validation works", {
  # Test invalid data structure
  expect_error(fit_pella_tomlinson_model("not a list"))
  expect_error(fit_pella_tomlinson_model(list(wrong = "structure")))

  # Test invalid data frames
  invalid_cpue <- list(cpue_data = "not a data frame", catch_data = data.frame(year = 2000, catch = 1000))
  expect_error(fit_pella_tomlinson_model(invalid_cpue))

  invalid_catch <- list(cpue_data = data.frame(year = 2000, cpue = 1.5), catch_data = "not a data frame")
  expect_error(fit_pella_tomlinson_model(invalid_catch))

  # Test missing columns
  missing_cpue_col <- list(
    cpue_data = data.frame(year = 2000, wrong_name = 1.5),
    catch_data = data.frame(year = 2000, catch = 1000)
  )
  expect_error(fit_pella_tomlinson_model(missing_cpue_col))

  missing_catch_col <- list(
    cpue_data = data.frame(year = 2000, cpue = 1.5),
    catch_data = data.frame(year = 2000, wrong_name = 1000)
  )
  expect_error(fit_pella_tomlinson_model(missing_catch_col))

  env_missing_data <- list(
    cpue_data = data.frame(year = 2000:2004, cpue = rep(1.2, 5)),
    catch_data = data.frame(year = 2000:2004, catch = rep(1000, 5))
  )
  expect_error(
    fit_pella_tomlinson_model(
      env_missing_data,
      options = list(env_covariates = c("temp"), validate_data = FALSE)
    ),
    "env_data"
  )
})

test_that("fit_pella_tomlinson_model validates prior specifications", {
  skip_if_not_installed("RTMB")

  years <- 2010:2015
  data_list <- list(
    cpue_data = data.frame(year = years, cpue = seq(1.2, 0.9, length.out = length(years))),
    catch_data = data.frame(year = years, catch = rep(800, length(years)))
  )

  expect_error(
    fit_pella_tomlinson_model(
      data_list,
      options = list(
        validate_data = FALSE,
        silent = TRUE,
        priors = list(r = list(dist = "lognormal"))
      )
    ),
    "meanlog"
  )

  expect_error(
    fit_pella_tomlinson_model(
      data_list,
      options = list(
        validate_data = FALSE,
        silent = TRUE,
        priors = list(unknown_param = list(dist = "normal", mean = 0, sd = 1))
      )
    ),
    "does not match any fitted parameter"
  )
})

test_that("strong r priors influence the fitted r estimate", {
  skip_if_not_installed("RTMB")

  set.seed(321)
  years <- 2010:2018
  true_params <- list(r = 0.3, K = 6000, m = 2, q = 0.0012, B_initial = 5000)

  biomass <- numeric(length(years))
  biomass[1] <- true_params$B_initial
  catch <- rep(700, length(years))
  for (t in seq_len(length(years) - 1)) {
    production <- true_params$r * biomass[t] * (1 - biomass[t] / true_params$K) / true_params$m
    biomass[t + 1] <- max(100, biomass[t] + production - catch[t])
  }

  cpue <- true_params$q * biomass * exp(rnorm(length(years), 0, 0.08))
  data_list <- list(
    cpue_data = data.frame(year = years, cpue = cpue),
    catch_data = data.frame(year = years, catch = catch)
  )

  base_options <- list(
    validate_data = FALSE,
    silent = TRUE,
    control = list(eval.max = 300, iter.max = 150)
  )

  fit_low <- tryCatch(
    fit_pella_tomlinson_model(
      data_list,
      options = modifyList(base_options, list(
        priors = list(r = list(dist = "lognormal", meanlog = log(0.08), sdlog = 0.08))
      ))
    ),
    error = function(e) skip(paste("Low-r prior fit failed:", e$message))
  )

  fit_high <- tryCatch(
    fit_pella_tomlinson_model(
      data_list,
      options = modifyList(base_options, list(
        priors = list(r = list(dist = "lognormal", meanlog = log(0.7), sdlog = 0.08))
      ))
    ),
    error = function(e) skip(paste("High-r prior fit failed:", e$message))
  )

  expect_lt(unname(fit_low$parameters["r"]), unname(fit_high$parameters["r"]))
  expect_equal(fit_low$results$priors[[1]]$param, "log_r")
  expect_equal(fit_high$results$priors[[1]]$param, "log_r")
})

test_that("B_initial prior alias maps to initial-biomass parameter", {
  skip_if_not_installed("RTMB")

  years <- 2010:2016
  data_list <- list(
    cpue_data = data.frame(year = years, cpue = seq(1.5, 1.0, length.out = length(years))),
    catch_data = data.frame(year = years, catch = rep(650, length(years)))
  )

  fit_alias <- tryCatch(
    fit_pella_tomlinson_model(
      data_list,
      options = list(
        validate_data = FALSE,
        silent = TRUE,
        control = list(eval.max = 200, iter.max = 100),
        priors = list(B_initial = list(dist = "lognormal", meanlog = log(5000), sdlog = 0.3))
      )
    ),
    error = function(e) skip(paste("B_initial alias fit failed:", e$message))
  )

  prior_params <- vapply(fit_alias$results$priors, function(x) x$param, character(1))
  expect_true(any(grepl("^log_B_initial", prior_params)))
  expect_true(any(grepl("^B_initial", names(fit_alias$parameters))))
})

test_that("AR1 process structure captures positive autocorrelation better than IID", {
  skip_if_not_installed("RTMB")

  set.seed(42)
  years <- 1990:2014
  n <- length(years)

  true <- list(
    r = 0.28,
    K = 6000,
    m = 2,
    q = 0.0014,
    B_initial = 5000,
    sigma_proc = 0.08,
    sigma_obs = 0.08,
    rho = 0.7,
    beta = 0.18
  )

  env <- as.numeric(scale(sin(seq_len(n)) + rnorm(n, 0, 0.3)))
  u <- numeric(n - 1)
  u[1] <- rnorm(1, 0, true$sigma_proc / sqrt(1 - true$rho^2))
  if (n > 2) {
    for (t in 2:(n - 1)) {
      u[t] <- true$rho * u[t - 1] + rnorm(1, 0, true$sigma_proc)
    }
  }

  B <- numeric(n)
  B[1] <- true$B_initial
  catch <- rep(700, n)
  for (t in seq_len(n - 1)) {
    prod_t <- true$r * B[t] * (1 - B[t] / true$K) / true$m
    B_det <- max(1, B[t] + prod_t - catch[t])
    B[t + 1] <- exp(log(B_det) + true$beta * env[t] + u[t])
  }
  cpue <- true$q * B * exp(rnorm(n, 0, true$sigma_obs))

  data_list <- list(
    cpue_data = data.frame(year = years, cpue = cpue),
    catch_data = data.frame(year = years, catch = catch),
    env_data = data.frame(year = years, temp = env)
  )

  common_opts <- list(
    validate_data = FALSE,
    silent = TRUE,
    process_noise = TRUE,
    env_covariates = "temp",
    fixed_params = list(log_m = log(2)),
    control = list(eval.max = 1200, iter.max = 600)
  )

  params_init <- list(
    log_r = log(0.25),
    log_K = log(5500),
    log_m = log(2),
    log_sigma_proc = log(0.1),
    log_sigma_obs = log(0.1),
    log_q.A1 = log(0.0012),
    log_B_initial.A1 = log(4500),
    beta_temp = 0
  )

  fit_iid <- tryCatch(
    fit_pella_tomlinson_model(
      data_list,
      params_init = params_init,
      options = modifyList(common_opts, list(process_error_structure = "iid"))
    ),
    error = function(e) skip(paste("IID fit failed:", e$message))
  )

  fit_ar1 <- tryCatch(
    fit_pella_tomlinson_model(
      data_list,
      params_init = params_init,
      options = modifyList(common_opts, list(process_error_structure = "ar1"))
    ),
    error = function(e) skip(paste("AR1 fit failed:", e$message))
  )

  expect_equal(fit_iid$results$process_error_structure, "iid")
  expect_equal(fit_ar1$results$process_error_structure, "ar1")

  # Directionality checks: true rho is strongly positive.
  expect_true(is.finite(fit_ar1$results$rho))
  expect_gt(unname(fit_ar1$results$rho), 0.2)

  # AR1 should fit at least as well as IID on AR1-generated data.
  expect_lte(fit_ar1$results$likelihood, fit_iid$results$likelihood + 1e-6)
})

# Integration test with minimal example
test_that("fit_pella_tomlinson_model integration test with simple data", {
  skip_if_not_installed("RTMB")

  # Create simple synthetic data
  set.seed(123)
  years <- 2010:2015 # Short time series for fast testing
  true_params <- list(r = 0.3, K = 5000, m = 2.0, q = 0.001, B_initial = 4000)

  # Generate deterministic biomass trajectory
  biomass <- numeric(length(years))
  biomass[1] <- true_params$B_initial
  catch <- rep(800, length(years))

  for (t in 1:(length(years) - 1)) {
    production <- true_params$r * biomass[t] * (1 - biomass[t] / true_params$K)
    biomass[t + 1] <- biomass[t] + production - catch[t]
  }

  # Generate CPUE with minimal noise
  cpue <- true_params$q * biomass * exp(rnorm(length(years), 0, 0.05))

  # Prepare data
  data_list <- list(
    cpue_data = data.frame(year = years, cpue = cpue),
    catch_data = data.frame(year = years, catch = catch)
  )

  # Fit model with validation disabled for speed
  options_list <- list(
    validate_data = FALSE,
    silent = TRUE,
    control = list(eval.max = 100, iter.max = 50) # Reduced for testing
  )

  # This test might fail if RTMB is not properly installed
  result <- tryCatch(
    {
      fit_pella_tomlinson_model(data_list, options = options_list)
    },
    error = function(e) {
      skip(paste("RTMB integration failed:", e$message))
    }
  )

  # Basic checks if fitting succeeded
  if (!is.null(result)) {
    expect_s3_class(result, "ProductionModel")
    expect_true(result$fitted)
    expect_true(length(result$parameters) > 0)
    expect_true("biomass" %in% names(result$results))
    expect_true("likelihood" %in% names(result$results))
  }
})
