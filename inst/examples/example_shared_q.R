# Two areas with a common CPUE catchability and a complete movement matrix.
library(SurplusProductionModel)
set.seed(42)
years <- 2000:2020
areas <- c("A1", "A2")
K <- c(A1 = 6000, A2 = 8000)
r <- 0.25
q <- 0.001
d0 <- 0.8
transition <- matrix(c(0.95, 0.05, 0.03, 0.97), nrow = 2, byrow = TRUE,
                     dimnames = list(from = areas, to = areas))
advance <- function(b, catch) {
  as.numeric((b + r * b * (1 - b / K) - catch) %*% transition)
}
equilibrium <- K
for (i in seq_len(50)) equilibrium <- advance(equilibrium, 0)
biomass <- matrix(NA_real_, length(years), 2, dimnames = list(years, areas))
biomass[1, ] <- d0 * equilibrium
catches <- cbind(A1 = seq(150, 350, length.out = length(years)),
                 A2 = seq(200, 450, length.out = length(years)))
for (i in seq_len(length(years) - 1)) {
  biomass[i + 1, ] <- advance(biomass[i, ], catches[i, ])
}
cpue_df <- do.call(rbind, lapply(seq_along(areas), function(a) {
  data.frame(year = years, area = areas[a],
             cpue = q * biomass[, a] * exp(rnorm(length(years), 0, 0.08)))
}))
catch_df <- do.call(rbind, lapply(seq_along(areas), function(a) {
  data.frame(year = years, area = areas[a], catch = catches[, a])
}))
# Do not include a label column, even if all observations have the same label.
processed <- prepare_model_data(cpue_df, catch_df)
starts <- prepare_starting_values(processed)
q_names <- grep("^log_q\\.", names(starts), value = TRUE)
starts$log_q_shared <- mean(unlist(starts[q_names]))
starts[q_names] <- NULL
fit_shared <- fit_pella_tomlinson_model(
  list(cpue_data = cpue_df, catch_data = catch_df,
       movement = list(transition_matrix = transition)),
  params_init = starts,
  options = list(shared_q = TRUE, spinup_years = 50,
                 fixed_params = list(log_m = log(2), log_d0 = log(d0)),
                 control = list(eval.max = 5000, iter.max = 2000))
)
print(fit_shared)
fit_shared$results$parameters["q_shared"]
plot_model_fit(fit_shared)
