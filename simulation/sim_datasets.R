
library(rstan)
library(dplyr)
library(bayesplot)
library(fields)
library(MASS)
library(raster)
library(sp)
library(sf)
library(geodist)

# ----------------------- Constants -------------------------- #

set.seed(123)

N_SIMS <- 50

OUTPUT_DIR <- 'sim_data'

AGG_FACTOR <- 16

# Range
Theta <- 200

# Marginal variance
Phi <- 1

# Outcome for locations
Beta_loc <- c(-0.25, 1.25, 1.5)

Beta_pos <- c(0.5, -1, 2)

# Additional spatial process
Theta_pos <- 100 # range
Phi_pos <- 0.5 # marginal variance

Beta_neg <-  c(3.5, 1.5, 2.5)

# Additional spatial process
Theta_neg <- 100 # range
Phi_neg <- 0.5 # marginal variance

# ------------- Covariates and Distance Matrix -------------- #

# add covariates from PRISM
caPr <- raster::stack("data/prism_pcas_ca.grd")
caPr.disc <- aggregate(caPr, fact = AGG_FACTOR)
plot(caPr.disc)

N <- length(caPr.disc[[1]][][!is.na(caPr.disc[[1]][])])
print(N)

# Covariates
X <- cbind(1, 
           caPr.disc[[1]][][!is.na(caPr.disc[[1]][])],
           caPr.disc[[2]][][!is.na(caPr.disc[[2]][])])

# get distance matrix
cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
coords <- xyFromCell(caPr.disc[[1]], cell = cells.all)
colnames(coords) <- c('longitude', 'latitude')

# Geodesic distances in kilometers
d <- geodist(coords, measure = 'geodesic' )/1000

# covariance function
cov_exp <- function(dist_matrix, alpha, rho){
  
  # Gaussian process covariance matrix
  cov_matrix <- matrix(NA, nrow=nrow(dist_matrix), ncol=ncol(dist_matrix))
  for (i in 1:nrow(dist_matrix)){
    for (j in 1:nrow(dist_matrix)){
      cov_matrix[i,j] <- (alpha ** 2) * exp(-dist_matrix[i,j] / rho)
    }
  }
  
  return(cov_matrix)
  
}

# Shared Latent Process Covariance
Sigma <- cov_exp(d, Phi, Theta)

print(max(Sigma[row(Sigma) == (col(Sigma) - 1)]))

# ----------------- Simulate Datasets ----------------------- #

summary_stats <- data.frame()

for (i in 1:N_SIMS){
  
  print(paste("iteration", i))
  
  # Keep generating datasets until prevalence and count criteria are met
  flag <- TRUE
  
  while (flag){
    
    # Shared latent process
    W <- mvrnorm(n = 1, mu = rep(0, N), Sigma)
    
    # ----------------- Simulate sampling effort ----------------- #
    
    # Randomly draw preferential sampling parameter
    Alpha_pos <- abs(rnorm(1, 0, 2))
    
    # Sampling counts
    Y_loc <- rpois(n = N, lambda = exp(X %*% Beta_loc + W))
    
    # ----------------- Disease positive counts ------------------ #
    
    # Indicators for whether each cell was observed at least once
    obs_inds <- Y_loc > 0
    
    # Additional spatial process
    Sigma_pos <- cov_exp(d, Phi_pos, Theta_pos)
    Eta_pos <- mvrnorm(n = 1, mu = rep(0, nrow(d)), Sigma_pos)
    
    # Disease positive counts over the entire study region
    Y_pos_all <- rpois(n = N, lambda = exp(X %*% Beta_pos + Eta_pos + Alpha_pos * W))
    
    # Observed disease positive counts
    Y_pos <- Y_pos_all[obs_inds]
    
    # ----------------- Disease negative counts ------------------ #
    
    # Preferential sampling parameter for disease negative counts
    Alpha_neg <- abs(rnorm(1, 0, 2))
    
    # Additional spatial process
    Sigma_neg <- cov_exp(d, Phi_neg, Theta_neg)
    Eta_neg <- mvrnorm(n = 1, mu = rep(0, nrow(d)), Sigma_neg)
    
    # Disease negative counts over the entire study region
    Y_neg_all <- rpois(n = N, lambda = exp(X %*% Beta_neg + Eta_neg + Alpha_neg * W))
    
    # Observed disease negative counts
    Y_neg <- Y_neg_all[obs_inds]
    
    # Check to see if criteria are met
    prev_i <- sum(Y_pos_all)/sum(Y_pos_all + Y_neg_all)
    frac_i <- mean(obs_inds)
    
    if (
      (prev_i > 0.02)
      & (prev_i < 0.4)
      & (sum(Y_pos) > 500)
      & (sum(Y_pos) < 50000)
      & (sum(Y_neg) > 500)
      & (sum(Y_neg) < 500000)
      & (frac_i > 0.1)
      & (frac_i < 0.7)
      & (Alpha_neg >= 0)
      & (Alpha_pos >= 0)
    ) {
      # Proceed to next dataset
      flag <- FALSE
    }
    
  }
  
  # ----------------------- Summary Statistics ------------------- #
  
  summary_stats <- rbind(
    summary_stats,
    data.frame(
      prevalence = sum(Y_pos)/(sum(Y_pos + Y_neg)),
      fraction_observed = round(mean(obs_inds), 3),
      num_pos_observed = sum(Y_pos),
      num_neg_observed = sum(Y_neg),
      alpha_p = Alpha_pos,
      alpha_n = Alpha_neg
    )
  )
  
  # ----------------------- Save Outputs ------------------------- #
  
  dat_i <- list(
    N=N,
    N_obs=sum(obs_inds),
    I_obs=(1:N)[obs_inds],
    P=3,
    X=X,
    D=d,
    Y_pos=Y_pos,
    Y_loc=Y_loc,
    Y_neg=Y_neg
  )
  
  params_i <- list(
    Theta = Theta,
    Phi = Phi,
    Beta_loc = Beta_loc,
    Beta_pos = Beta_pos,
    Theta_pos = Theta_pos,
    Phi_pos = Phi_pos,
    Beta_neg = Beta_neg,
    Theta_neg = Theta_neg,
    Phi_neg = Phi_neg,
    Alpha_pos = Alpha_pos,
    Alpha_neg = Alpha_neg,
    W = W,
    Eta_pos = Eta_pos,
    Eta_neg = Eta_neg
  )
  
  counts_i <- list(
    Y_pos_all = Y_pos_all,
    Y_neg_all = Y_neg_all
  )
  
  saveRDS(dat_i, file = paste0(OUTPUT_DIR, '/dataset_', i, '.rds'))
  saveRDS(params_i, file = paste0(OUTPUT_DIR, '/parameters_', i, '.rds'))
  saveRDS(counts_i, file = paste0(OUTPUT_DIR, '/counts_', i, '.rds'))
  
}

# -------------------- Check Data Characteristics -------------------- #

# Prevalence
hist(summary_stats$prevalence)
summary(summary_stats$prevalence)

# Numbers observed
summary(summary_stats$num_pos_observed)
summary(summary_stats$num_neg_observed)

# Difference in preferential sampling params
hist(summary_stats$alpha_p - summary_stats$alpha_n)

# Fractions observed
hist(summary_stats$fraction_observed)
summary(summary_stats$fraction_observed)

write.csv(summary_stats, 'sim_summary_stats.csv')
