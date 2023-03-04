
library(rstan)
library(sp)
library(raster)
library(MASS)
library(posterior)
library(dplyr)

set.seed(123)

# -------------------- Load ------------------------------------- #

dat <- readRDS(file = 'data/plague_data_agg7.rds')
mod_fit <- readRDS(file = 'cdph_fits/plague_slp_fit_10.rds')

caPr <- raster::stack("data/prism_pcas_ca.grd")
caPr.disc <- aggregate(caPr, fact = 7)
plot(caPr.disc)

N <- length(caPr.disc[[1]][][!is.na(caPr.disc[[1]][])])
print(N)

# -------------------- Simulate Data ---------------------------- #

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
offset <- dat$Y_loc[dat$I_obs]
exp_pos <- offset * exp(X[dat$I_obs,] %*% beta_pos_hat + eta_pos_hat + alpha_pos_hat * w_hat[dat$I_obs])
exp_neg <- offset * exp(X[dat$I_obs,] %*% beta_neg_hat + eta_neg_hat + alpha_neg_hat * w_hat[dat$I_obs])

# Simulate values
Y_loc <- dat$Y_loc
Y_pos_sim <- rpois(n = length(dat$I_obs), lambda = exp_pos)
Y_neg_sim <- rpois(n = length(dat$I_obs), lambda = exp_neg)

# -------------------- Inspect distributions -------------------- #
print(length(Y_pos_sim))
print(length(Y_neg_sim))
print(length(Y_loc))

hist(Y_pos_sim)
hist(Y_neg_sim)
hist(Y_loc)

print(summary(Y_pos_sim))
print(summary(Y_neg_sim))
print(sum(Y_pos_sim)/sum(Y_pos_sim + Y_neg_sim))

# --------------------------- Save ------------------------------ #

dat_sim <- list(
  Y_pos = Y_pos_sim,
  Y_neg = Y_neg_sim,
  Y_loc = Y_loc,
  X = dat$X,
  prism_rasters = caPr.disc
)

saveRDS(dat_sim, 'plague_data_simulated.rds')

# -------------------- Create CSV Version ----------------------- #

# all raster cells that have values
cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
cells_obs <- cells.all[dat$I_obs]

# cell ids and counts
tmp <- data.frame(
  cell_id = cells_obs,
  Y_pos = Y_pos_sim,
  Y_neg = Y_neg_sim)

# all ids and sampling counts
sim_df <- data.frame(
  cell_id = cells.all,
  Y_loc = Y_loc)

# combine
sim_df <- sim_df %>%
  left_join(tmp,
            by = 'cell_id') %>%
  mutate(Y_pos = ifelse(is.na(Y_pos), 0, Y_pos),
         Y_neg = ifelse(is.na(Y_neg), 0, Y_neg))

# prism variables
sim_df$X_pc1 = X[,2]
sim_df$X_pc2 = X[,3]

write.csv(sim_df, 'sim_df.csv', row.names = F)

# -------------------- Document --------------------------------- #

prism_rasters = dat_sim$prism_rasters
Y_loc = dat_sim$Y_loc
Y_pos = dat_sim$Y_pos

# use the first raster layer for the study region
study_region <- prism_rasters[[1]]

# overlay counts of sampling events
sampling_raster <- study_region
sampling_raster[][!is.na(sampling_raster[])] <- Y_loc
plot(sampling_raster)

# the counts of Y_pos are indexed with respect to cells that have been sampled by the
# surveillance system, not all cells, so we have to apply additional logic
rodent_raster <- study_region
is_sampled <- Y_loc > 0
is_in_region <- !is.na(rodent_raster[])
# fill in the unobserved cells with 0
rodent_raster[][!is.na(rodent_raster[])] <- 0
rodent_raster[][is_in_region][is_sampled] <- Y_pos
plot(rodent_raster)
