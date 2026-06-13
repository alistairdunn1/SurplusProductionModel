# Consistency between the standard Pella-Tomlinson production dynamics and the
# canonical reference-point formulae. These tests guard against the production
# normalisation and the reference points drifting apart again.

test_that("reference points match the numerical optimum of the production function", {
  for (m in c(1.5, 2, 3, 4)) {
    r <- 0.3
    K <- 5000
    rp <- pella_tomlinson_reference_points(r = r, K = K, m = m)

    # Numerically maximise production over biomass
    opt <- optimize(
      function(B) SurplusProductionModel:::.pt_production(B, r, K, m),
      interval = c(1e-6 * K, K), maximum = TRUE
    )

    expect_equal(rp$bmsy, opt$maximum, tolerance = 1e-4)
    expect_equal(rp$msy, opt$objective, tolerance = 1e-4)
    expect_equal(rp$fmsy, rp$msy / rp$bmsy, tolerance = 1e-10)
  }
})

test_that("Schaefer (m = 2) reduces to the textbook reference points", {
  rp <- pella_tomlinson_reference_points(r = 0.3, K = 5000, m = 2)
  expect_equal(rp$bmsy, 2500)
  expect_equal(rp$fmsy, 0.15)
  expect_equal(rp$msy, 375)
})

test_that("Fox (m = 1) limit gives K/e, r, rK/e", {
  rp <- pella_tomlinson_reference_points(r = 0.3, K = 5000, m = 1)
  expect_equal(rp$bmsy, 5000 / exp(1), tolerance = 1e-10)
  expect_equal(rp$fmsy, 0.3, tolerance = 1e-10)
  expect_equal(rp$msy, 0.3 * 5000 / exp(1), tolerance = 1e-10)
})

test_that("equilibrium harvest rate at BMSY equals FMSY", {
  for (m in c(1, 1.5, 2, 3)) {
    r <- 0.25
    K <- 8000
    rp <- pella_tomlinson_reference_points(r = r, K = K, m = m)
    f_at_bmsy <- pt_equilibrium_f(r = r, K = K, m = m, biomass = rp$bmsy)
    expect_equal(f_at_bmsy, rp$fmsy, tolerance = 1e-8)
  }
})

test_that("production at BMSY equals MSY", {
  for (m in c(1.5, 2, 3)) {
    r <- 0.4
    K <- 1000
    rp <- pella_tomlinson_reference_points(r = r, K = K, m = m)
    prod_at_bmsy <- SurplusProductionModel:::.pt_production(rp$bmsy, r, K, m)
    expect_equal(prod_at_bmsy, rp$msy, tolerance = 1e-6)
  }
})
