# Quick test: does tmbstan work with RTMB objects?
library(RTMB)
library(tmbstan)

x <- c(1.2, 0.8, 1.0)
obj <- MakeADFun(function(p) {
  -sum(dnorm(x, p$mu, exp(p$log_sd), log = TRUE))
}, list(mu = 0, log_sd = 0), silent = TRUE)

cat("obj class:", class(obj), "\n")
cat("obj$par:", obj$par, "\n")
cat("names:", names(obj$par), "\n")

# Try tmbstan with very short run
fit <- tmbstan(obj, chains = 1, iter = 200, warmup = 100, seed = 42)
cat("stanfit class:", class(fit), "\n")
cat("summary:\n")
print(summary(fit)$summary)
