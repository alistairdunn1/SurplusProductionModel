# Test RTMB 'map' argument for fixing parameters
library(RTMB)

set.seed(42)
n <- 20
x_true <- cumsum(rnorm(n, 0, 0.3))
y_obs <- x_true + rnorm(n, 0, 0.5)
y <- y_obs

# Same model
f <- function(parms) {
  getAll(parms, warn = FALSE)
  nll <- 0
  for (t in 2:length(x)) {
    nll <- nll - dnorm(x[t], x[t - 1], exp(log_sd_proc), log = TRUE)
  }
  for (t in 1:length(y)) {
    nll <- nll - dnorm(y[t], x[t], exp(log_sd_obs), log = TRUE)
  }
  nll
}

parms <- list(
  x = rep(0, n),
  log_sd_proc = log(0.3),
  log_sd_obs = log(0.5)
)

# Fix log_sd_proc using map (factor with NA = fixed)
obj <- MakeADFun(f, parms,
  random = "x",
  map = list(log_sd_proc = factor(NA)),
  silent = TRUE
)
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep <- sdreport(obj)

cat("Fixed params (with log_sd_proc fixed):", opt$par, "\n")
cat("Convergence:", opt$convergence, "\n")
summ_fixed <- summary(rep, "fixed")
cat("\nFixed param estimates with SE:\n")
print(summ_fixed)

# Test ADREPORT for derived quantities
f2 <- function(parms) {
  getAll(parms, warn = FALSE)
  nll <- 0
  sd_proc <- exp(log_sd_proc)
  sd_obs <- exp(log_sd_obs)
  for (t in 2:length(x)) {
    nll <- nll - dnorm(x[t], x[t - 1], sd_proc, log = TRUE)
  }
  for (t in 1:length(y)) {
    nll <- nll - dnorm(y[t], x[t], sd_obs, log = TRUE)
  }
  # Report derived quantities
  ADREPORT(sd_proc)
  ADREPORT(sd_obs)
  nll
}

obj2 <- MakeADFun(f2, parms, random = "x", silent = TRUE)
opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr)
rep2 <- sdreport(obj2)

cat("\nADREPORT derived quantities:\n")
summ_report <- summary(rep2, "report")
print(summ_report)
