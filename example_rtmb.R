


library(RTMB)
library(TMBhelper)
library(tmbstan)
library(SparseNuts)

# number of indices
n_q <- 1L

# time
n_time <- 30

# parameters and initial values
parameters <- list(

    log_r        = log(0.35), # Growth rate
    log_K        = log(1500), # Carrying capacity
    log_z        = log(1), # With z=1, the PT model is the Shaeffer model
    log_q        = log(rep(0.001, n_q)), # catchability
    log_cpue_pow = log(1),
    log_cpue_pro = log(rep(1e-6, n_q)), # process error for CPUE
    log_sigmap   = log(0.1), # process error
    log_B        = log(rep(1500, n_time - 1))
)

priors <- list()
priors[["log_r"]]        <- list(type = "normal", par1 = log(0.35), par2 = 0.5, index = which("log_r"        == names(parameters)))
priors[["log_K"]]        <- list(type = "normal", par1 = log(1000), par2 = 1,   index = which("log_K"        == names(parameters)))
priors[["log_z"]]        <- list(type = "normal", par1 = 0,         par2 = 0.1, index = which("log_z"        == names(parameters)))
priors[["log_q"]]        <- list(type = "normal", par1 = -5,        par2 = 2,   index = which("log_q"        == names(parameters)))
priors[["log_cpue_pow"]] <- list(type = "normal", par1 = 0,         par2 = 0.5, index = which("log_cpue_pow" == names(parameters)))
priors[["log_cpue_pro"]] <- list(type = "normal", par1 = 0,         par2 = 0.5, index = which("log_cpue_pro" == names(parameters)))
priors[["log_sigmap"]]   <- list(type = "normal", par1 = 0,         par2 = 0.5, index = which("log_sigmap"   == names(parameters)))

priors

# Model
fun <- function(parameters, data) {
    
    # load all available parameters
    # and data
    getAll(parameters, data, warn = FALSE)
    
    # back transform estimated
    # values
    r <- exp(log_r)
    K <- exp(log_K)
    z <- exp(log_z)
    B <- c(K, exp(log_B))
    
    n_time <- length(catch_time)
    
    production <- r / z * B * (1 - (B / K)^z)
    
    Bpred <- numeric(n_time)
    Bpred[1] <- K
    for (t in 2:n_time) {
        Bpred[t] <- B[t - 1] + production[t - 1] - catch_obs[t - 1]
    }
    
    nll_B <- -sum(dnorm(log(B), mean = log(Bpred), sd = exp(log_sigmap), log = TRUE))
    REPORT(B)
    
    cpue_obs <- OBS(cpue_obs)
    cpue_pow <- exp(log_cpue_pow)
    
    cpue_pred <- exp(log_q) * B[cpue_time]
    nll_cpue  <- -dlnorm(x = cpue_obs, meanlog = log(cpue_pred), sdlog = log_sigmao, log = TRUE)
    REPORT(cpue_pred)
    
    nll <- nll_B + sum(nll_cpue)
    return(nll)
}

Lwr <- rep(-Inf, length(obj$par))
Upr <- rep(Inf, length(obj$par))
Lwr[grep("log_r", names(obj$par))] <- -10
Upr[grep("log_r", names(obj$par))] <- 10
Lwr[grep("log_K", names(obj$par))] <- log(300)
Upr[grep("log_K", names(obj$par))] <- log(150000)
Lwr[grep("log_z", names(obj$par))] <- log(0)
Upr[grep("log_z", names(obj$par))] <- log(2)
Lwr[grep("log_q", names(obj$par))] <- -25
Upr[grep("log_q", names(obj$par))] <- 1
Lwr[grep("log_sigmap", names(obj$par))] <- log(0)
Upr[grep("log_sigmap", names(obj$par))] <- log(2)

map <- list(
    log_z = factor(NA)
)

cmb <- function(f, d) function(p) f(p, d)

obj <- MakeADFun(cmb(fun, data), parameters, random = "log_B", map = map)



control <- list(eval.max = 100000, iter.max = 100000)

opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr, control = control, lower = Lwr, upper = Upr)
opt <- nlminb(start = opt$par, objective = obj$fn, gradient = obj$gr, control = control, lower = Lwr, upper = Upr)
# check_estimability(obj = obj)


