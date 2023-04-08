library(rstan)
library(paws)

# ----------------------- Constants ------------------------- #

OUT_FILE <- 'plague_spatial_poisson_fit_spois.rds'

DATA_FILE <- 'plague_data_agg7.rds'

set.seed(123)

Sys.setenv(
  AWS_ACCESS_KEY_ID = "",
  AWS_SECRET_ACCESS_KEY = "",
  AWS_REGION = ""
)

svc <- paws::s3()

# ----------------------- Analysis ------------------------- #

# retrieve data
svc$download_file(
  Bucket = "plague-analysis",
  Key = DATA_FILE,
  Filename = DATA_FILE
)
dat <- readRDS(file = DATA_FILE)

# Subset data for Spatial Poisson model
dat_spois <- list(
  N_obs = dat$N_obs,
  P = 3,
  X = dat$X[dat$I_obs,],
  D = dat$D[dat$I_obs, dat$I_obs],
  Y_pos = dat$Y_pos,
  Y_neg = dat$Y_neg,
  N_samp = dat$Y_loc[dat$I_obs]
)

print(paste("fraction observed:", round(dat$N_obs/dat$N, 2)))
print(paste("Y (+) total:", sum(dat$Y_pos)))
print(paste("Y (-) total:", sum(dat$Y_neg)))
print(paste("N: ", dat$N))

# Fit STAN model
fit <- stan(
  file = 'spatial_poisson.stan', 
  data = dat_spois,
  chains = 4,
  cores = 4)

saveRDS(fit, file = OUT_FILE)

svc$put_object(
  Bucket = "plague-analysis",
  Body = OUT_FILE,
  Key = OUT_FILE
)
