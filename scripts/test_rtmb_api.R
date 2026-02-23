# Quick test of RTMB API for state-space model pattern
library(RTMB)

set.seed(42)
n <- 20
x_true <- cumsum(rnorm(n, 0, 0.3))
y_obs <- x_true + rnorm(n, 0, 0.5)

# RTMB way: define function taking 'parameters' list
f <- function(parms) {
  getAll(parms, warn = FALSE)
  nll <- 0
  # Random effects (states)
  for (t in 2:length(x)) {
    nll <- nll - dnorm(x[t], x[t - 1], exp(log_sd_proc), log = TRUE)
  }
  # Observations
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

# Fixed data (must exist in environment)
y <- y_obs

obj <- MakeADFun(f, parms, random = "x", silent = TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep <- sdreport(obj)

cat("Fixed params:", opt$par, "\n")
cat("Convergence:", opt$convergence, "\n")
cat("sdreport class:", class(rep), "\n")

# Extract random effect estimates
summ <- summary(rep, "random")
cat("\nFirst 3 state estimates:\n")
print(head(summ, 3))

# Extract fixed parameter SEs
summ_fixed <- summary(rep, "fixed")
cat("\nFixed param estimates with SE:\n")
print(summ_fixed)
