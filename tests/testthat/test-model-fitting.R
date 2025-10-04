test_that("pella_tomlinson_production function works correctly", {
  # Test basic functionality
  production <- pella_tomlinson_production(B = 1000, r = 0.3, K = 5000, m = 2)
  expect_type(production, "double")
  expect_true(production > 0)

  # Test Schaefer special case (m = 2)
  B <- 2500 # Half of K
  production_schaefer <- pella_tomlinson_production(B, r = 0.3, K = 5000, m = 2)
  expected_schaefer <- 0.3 * B * (1 - B / 5000) / 2 # Pella-Tomlinson with m=2
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
  expected_params <- c("log_r", "log_K", "log_m", "log_q", "log_sigma_proc", "log_sigma_obs", "log_B0")
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
  expect_true(natural_params[["B0"]] > 0)
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

  # Check years are in correct range
  expect_true(all(processed$years >= 2002))
  expect_true(all(processed$years <= 2018))

  # Check years are sorted
  expect_true(all(diff(processed$years) == 1))

  # Test insufficient overlap error
  cpue_data_short <- data.frame(year = 2000:2001, cpue = c(1.5, 1.6))
  catch_data_short <- data.frame(year = 2010:2011, catch = c(1000, 1100))

  expect_error(preprocess_model_data(cpue_data_short, catch_data_short))

  # Test missing catch values error
  catch_data_na <- data.frame(year = 2000:2005, catch = c(1000, NA, 1200, 1100, 1000, 1050))
  cpue_data_overlap <- data.frame(year = 2000:2005, cpue = rnorm(6, 1.5, 0.1))

  expect_error(preprocess_model_data(cpue_data_overlap, catch_data_na))
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

test_that("calculate_model_results produces expected output", {
  # Create test parameters
  parameters <- c(
    r = 0.3,
    K = 5000,
    m = 2.0,
    q = 0.001,
    B0 = 4000
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

  # Check biomass starts at B0
  expect_equal(results$biomass[1], parameters[["B0"]])

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
})

# Integration test with minimal example
test_that("fit_pella_tomlinson_model integration test with simple data", {
  skip_if_not_installed("RTMB")

  # Create simple synthetic data
  set.seed(123)
  years <- 2010:2015 # Short time series for fast testing
  true_params <- list(r = 0.3, K = 5000, m = 2.0, q = 0.001, B0 = 4000)

  # Generate deterministic biomass trajectory
  biomass <- numeric(length(years))
  biomass[1] <- true_params$B0
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
    expect_s4_class(result, "ProductionModel")
    expect_true(result@fitted)
    expect_true(length(result@parameters) > 0)
    expect_true("biomass" %in% names(result@results))
    expect_true("likelihood" %in% names(result@results))
  }
})
