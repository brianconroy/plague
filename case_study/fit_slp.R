library(rstan)
library(paws)

# ----------------------- Constants ------------------------- #

OUT_FILE <- 'cdph_fits/plague_slp_fit_agg7.rds'

DATA_FILE <- 'data/plague_data_agg7.rds'

set.seed(123)

# Fill these in
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

print(paste("fraction observed:", round(dat$N_obs/dat$N, 2)))
print(paste("Y (+) total:", sum(dat$Y_pos)))
print(paste("Y (-) total:", sum(dat$Y_neg)))
print(paste("N:", sum(dat$N)))

# Fit STAN model
fit <- stan(
  file = 'slp.stan',
  data = dat,
  chains = 2,
  cores = 2)

saveRDS(fit, file = OUT_FILE)

svc$put_object(
  Bucket = "plague-analysis",
  Body = OUT_FILE,
  Key = OUT_FILE
)
