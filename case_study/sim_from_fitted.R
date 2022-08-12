
library(rstan)
library(sp)
library(raster)
library(MASS)
library(posterior)

set.seed(123)

dat <- readRDS(file = 'data/plague_data_agg7.rds')
mod_fit <- readRDS(file = 'cdph_fits/plague_slp_fit_agg7.rds')

caPr <- raster::stack("data/prism_pcas_ca.grd")
caPr.disc <- aggregate(caPr, fact = 7)
plot(caPr.disc)

N <- length(caPr.disc[[1]][][!is.na(caPr.disc[[1]][])])
print(N)

# Get posterior samples
posterior <- as.matrix(mod_fit)
beta_pos_samples <- posterior[,grepl('beta_pos', colnames(posterior))]
beta_neg_samples <- posterior[,grepl('beta_neg', colnames(posterior))]
beta_loc_samples <- posterior[,grepl('beta_loc', colnames(posterior))]
alpha_pos_samples <- posterior[,'alpha_pos']
alpha_neg_samples <- posterior[,'alpha_neg']
print(paste("N col w:", ncol(posterior[,grepl('w', colnames(posterior))])))
w_samples <- posterior[,grepl('w', colnames(posterior))]
eta_pos_samples <- posterior[,
                             grepl('eta_pos', colnames(posterior))
                             & !(grepl('beta', colnames(posterior)))
                             & !(grepl('theta', colnames(posterior)))]
eta_neg_samples <- posterior[,
                             grepl('eta_neg', colnames(posterior))
                             & !(grepl('beta', colnames(posterior)))
                             & !(grepl('theta', colnames(posterior)))]

# Posterior means
eta_pos_hat <- colMeans(eta_pos_samples)
eta_neg_hat <- colMeans(eta_neg_samples)
alpha_pos_hat <- mean(alpha_pos_samples)
alpha_neg_hat <- mean(alpha_neg_samples)
beta_pos_hat <- colMeans(beta_pos_samples)
beta_neg_hat <- colMeans(beta_neg_samples)
beta_loc_hat <- colMeans(beta_loc_samples)
w_hat <- colMeans(w_samples)

# Calculate expected values of disease positive and negative specimen, and sampling effort
X <- dat$X
exp_pos <- exp(X[dat$I_obs,] %*% beta_pos_hat + eta_pos_hat + alpha_pos_hat * w_hat[dat$I_obs])
exp_neg <- exp(X[dat$I_obs,] %*% beta_neg_hat + eta_neg_hat + alpha_neg_hat * w_hat[dat$I_obs])

# Simulate values
Y_loc <- dat$Y_loc
Y_pos_sim <- rpois(n = length(dat$I_obs), lambda = exp_pos)
Y_neg_sim <- rpois(n = length(dat$I_obs), lambda = exp_neg)

print(length(Y_pos_sim))
print(length(Y_neg_sim))
print(length(Y_loc))

hist(Y_pos_sim)
hist(Y_neg_sim)
hist(Y_loc)

dat_sim <- list(
  Y_pos = Y_pos_sim,
  Y_neg = Y_neg_sim,
  Y_loc = Y_loc,
  X = dat$X,
  prism_rasters = caPr.disc
)

saveRDS(dat_sim, 'plague_data_simulated.rds')
