


library(RTMB)
library(tmbstan)

# time
n_time <- 30

# parameters and initial values
parameters <- list(

    log_r        = log(0.35), # Growth rate
    log_K        = log(1500), # Carrying capacity
    log_q        = log(0.001), # catchability
    log_B        = log(rep(1500, n_time - 1))
)

input_data        <- list()
input_data$year   <- 1:n_time
input_data$catch  <- seq(300, 400, length = n_time) 
input_data$index  <- seq(2, 1, length = n_time)
input_data$sigmao <- 0.1
input_data$sigmap <- 0.01
input_data$p      <- 1

# Model
fun <- function(parameters, data) {
    
    # load all available parameters
    # and data
    getAll(parameters, data, warn = FALSE)
    
    # back transform estimated
    # values
    r       <- exp(log_r)
    K       <- exp(log_K)
    B       <- c(K, exp(log_B))
    
    # time dimension
    n_time <- length(year)
    
    # observed input values
    index  <- OBS(index)
    #sigmao <- OBS(sigmao)
    #sigmap <- OBS(sigmap)
    
    # error terms
    sigmap2 <- sigmap^2
    sigmao2 <- sigmao^2
    
    # AD vectors
    production <- AD(numeric(n_time))
    mu         <- AD(numeric(n_time))
    
    # surplus productive
    production <- r / p * B * (1 - (B / K)^p)

    # expected annual dynamics
    mu[1] <- B[1]
    for (t in 2:n_time) {
        mu[t] <- B[t - 1] + production[t - 1] - catch[t - 1]
    }
    
    # state equation
    B %~% dlnorm(log(mu) - sigmap2 / 2, sigmap)
    
    # expected index
    index_pred <- exp(log_q) * B
    
    # observation equation
    index_pred %~% dlnorm(log(index) - sigmao2 / 2, sigmao)
    
    # prior
    log_r %~% dnorm(log(0.35), 0.1)
    
    # AD reports
    REPORT(B)
    REPORT(index_pred)
}

# set up objective function
cmb <- function(f, d) function(p) f(p, d)
obj <- MakeADFun(cmb(fun, input_data), parameters, random = "log_B")

# parameter bounds
Lwr <- rep(-Inf, length(obj$par))
Upr <- rep(Inf, length(obj$par))
Lwr[grep("log_K", names(obj$par))] <- log(300)
Upr[grep("log_K", names(obj$par))] <- log(150000)
Lwr[grep("log_q", names(obj$par))] <- -25
Upr[grep("log_q", names(obj$par))] <- 1

# fit
control <- list(eval.max = 100000, iter.max = 100000)

opt <- nlminb(start = obj$par, objective = obj$fn, gradient = obj$gr, control = control, lower = Lwr, upper = Upr)
opt <- nlminb(start = opt$par, objective = obj$fn, gradient = obj$gr, control = control, lower = Lwr, upper = Upr)

# report
summary(sdreport(obj))

# run stan
opt_stan <- tmbstan(obj, init = function() parameters)
traceplot(opt_stan, pars = names(opt$par))


