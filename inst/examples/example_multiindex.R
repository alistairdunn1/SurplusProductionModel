## Example: Multi-area with multiple CPUE indices per area
set.seed(1)

years <- 2000:2010
areas <- c("A1", "A2")
labels <- c("LL", "Trawl")
nY <- length(years)
nA <- length(areas)
nL <- length(labels)

r <- 0.25
K <- 6000
m <- 2
B_initial <- c(A1 = 3000, A2 = 3200)
q_mat <- matrix(
  c(
    0.0004, 0.00025,
    0.0005, 0.0003
  ),
  nrow = 2, byrow = TRUE,
  dimnames = list(area = areas, label = labels)
)

B <- matrix(NA_real_, nY, nA, dimnames = list(year = years, area = areas))
B[1, ] <- B_initial
C <- matrix(0, nY, nA, dimnames = list(year = years, area = areas))
C[, "A1"] <- 180
C[, "A2"] <- 220
for (t in 1:(nY - 1)) {
  P <- r * B[t, ] * (1 - (B[t, ] / K))
  B[t + 1, ] <- pmax(B[t, ] + P - C[t, ], 0.01)
}

sigma_obs <- 0.25
cpue_arr <- array(NA_real_, dim = c(nY, nA, nL), dimnames = list(year = years, area = areas, label = labels))
for (a in 1:nA) {
  for (l in 1:nL) {
    mu <- q_mat[a, l] * B[, a]
    cpue_arr[, a, l] <- mu * exp(rnorm(nY, 0, sigma_obs))
  }
}

cpue_df <- do.call(rbind, lapply(1:nA, function(a) {
  do.call(rbind, lapply(1:nL, function(l) {
    data.frame(year = years, area = areas[a], label = labels[l], cpue = cpue_arr[, a, l])
  }))
}))
catch_df <- do.call(rbind, lapply(1:nA, function(a) {
  data.frame(year = years, area = areas[a], catch = C[, a])
}))

data_list <- list(cpue_data = cpue_df, catch_data = catch_df)

fit <- SurplusProductionModel::fit_pella_tomlinson_model(
  data_list,
  options = list(
    fixed_params = list(log_m = log(m)),
    control = list(eval.max = 5000, iter.max = 2000)
  )
)
print(fit)
plots <- SurplusProductionModel::plot_model_fit(fit)
