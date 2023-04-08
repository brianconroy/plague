library(rstan)

# ----------------------- Constants ------------------------- #

set.seed(123)

OUT_FILE <- 'cdph_fits/plague_poisson_fit_pois.rds'

DATA_FILE <- 'data/plague_data_agg7.rds'

# ----------------------- Analysis ------------------------- #

dat <- readRDS(file = DATA_FILE)

# Subset data for Spatial Poisson model
dat_pois <- list(
  N_obs = dat$N_obs,
  P = 3,
  X = dat$X[dat$I_obs,],
  Y_pos = dat$Y_pos,
  Y_neg = dat$Y_neg,
  N_samp = dat$Y_loc[dat$I_obs]
)

print(paste("fraction observed:", round(dat$N_obs/dat$N, 2)))
print(paste("Y (+) total:", sum(dat$Y_pos)))
print(paste("Y (-) total:", sum(dat$Y_neg)))

# Fit STAN model
fit <- stan(
  file = 'poisson.stan', 
  data = dat_pois,
  chains = 4,
  cores = 4)

saveRDS(fit, file = OUT_FILE)
