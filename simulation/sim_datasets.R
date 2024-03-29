
# This script simulates disease surveillance datasets
# from the proposed shared latent process model.

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

AGG_FACTOR <- 15

# Range
Theta <- 100

# Marginal variance
Phi <- 2

# Outcome for locations
Beta_loc <- c(-1, 1.5, -1)

Beta_pos <- c(-1, 2.5, 1.5)

# Additional spatial process
Theta_pos <- 100 #  Range
Phi_pos <- 0.5   #  Marginal variance

Beta_neg <-  c(1, 0.5, 0.75)

# Additional spatial process
Theta_neg <- 100 #  Range
Phi_neg <- 0.5   #  Marginal variance

# ------------- Covariates and Distance Matrix -------------- #

# Load covariates from PRISM
caPr <- raster::stack("data/prism_pcas_ca.grd")
caPr.disc <- aggregate(caPr, fact = AGG_FACTOR)
plot(caPr.disc)

N <- length(caPr.disc[[1]][][!is.na(caPr.disc[[1]][])])
print(N)

# Covariates
X <- cbind(1, 
           caPr.disc[[1]][][!is.na(caPr.disc[[1]][])],
           caPr.disc[[2]][][!is.na(caPr.disc[[2]][])])

# Get distance matrix
cells.all <- c(1:ncell(caPr.disc))[!is.na(values(caPr.disc[[1]]))]
coords <- xyFromCell(caPr.disc[[1]], cell = cells.all)
colnames(coords) <- c('longitude', 'latitude')

# Geodesic distances in kilometers
d <- geodist(coords, measure = 'geodesic')/1000

# Covariance function
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
  
  counter <- 0
  
  while (flag){
    
    counter <- counter + 1
    
    # Shared latent process
    W <- mvrnorm(n = 1, mu = rep(0, N), Sigma)
    
    # ----------------- Simulate sampling effort ----------------- #
    
    # Randomly draw preferential sampling parameter
    Alpha_pos <- abs(rnorm(1, 0, 1))
    
    # Sampling counts
    Y_loc <- rpois(n = N, lambda = exp(X %*% Beta_loc + W))
    
    # ----------------- Disease positive counts ------------------ #
    
    # Indicators for whether each cell was observed at least once
    obs_inds <- Y_loc > 0
    
    # Additional spatial process
    Sigma_pos <- cov_exp(d, Phi_pos, Theta_pos)
    Eta_pos <- mvrnorm(n = 1, mu = rep(0, nrow(d)), Sigma_pos)
    
    # Offset
    N_offset <- Y_loc
    
    # Disease positive counts over the entire study region
    Y_pos_all <- rpois(n = N, lambda = N_offset * exp(X %*% Beta_pos + Eta_pos + Alpha_pos * W))
    
    # Observed disease positive counts
    Y_pos <- Y_pos_all[obs_inds]
    
    # ----------------- Disease negative counts ------------------ #
    
    # Preferential sampling parameter for disease negative counts
    Alpha_neg <- abs(rnorm(1, 0, 1))
    
    # Additional spatial process
    Sigma_neg <- cov_exp(d, Phi_neg, Theta_neg)
    Eta_neg <- mvrnorm(n = 1, mu = rep(0, nrow(d)), Sigma_neg)
    
    # Disease negative counts over the entire study region
    Y_neg_all <- rpois(n = N, lambda = N_offset * exp(X %*% Beta_neg + Eta_neg + Alpha_neg * W))
    
    # Observed disease negative counts
    Y_neg <- Y_neg_all[obs_inds]
    
    # Check to see if criteria are met
    prev_i <- sum(Y_pos_all)/sum(Y_pos_all + Y_neg_all)
    frac_i <- mean(obs_inds)

    if (
      (prev_i > 0.02)
      & (prev_i < 0.3)
      & (sum(Y_pos) > 500)
      & (sum(Y_pos) < 50000)
      & (sum(Y_neg) > 500)
      & (sum(Y_neg) < 500000)
      & (max(Y_pos) < 5000)
      & (max(Y_neg) < 10000)
      & (frac_i > 0.33)
      & (frac_i < 0.66)
      & (Alpha_neg <= 2)
      & (Alpha_pos <= 2)
    ) {
      
      # Proceed to next dataset
      flag <- FALSE
      
      par(mfrow = c(3,2))
      plot(x = abs(X %*% Beta_pos),
           y = abs(Alpha_pos * W),
           main = 'Y(+): Covariates vs Alpha*W')
      abline(0, 1)
      
      plot(x = abs(X %*% Beta_neg),
           y = abs(Alpha_neg * W),
           main = 'Y(-): Covariates vs Alpha*W')
      abline(0, 1)
      
      plot(x = abs(Eta_pos),
           y = abs(Alpha_pos * W),
           main = 'Y(+): Eta vs Alpha*W')
      abline(0, 1)
      
      plot(x = abs(Eta_neg),
           y = abs(Alpha_neg * W),
           main = 'Y(-): Eta vs Alpha*W')
      abline(0, 1)
      
      hist(Alpha_pos * W - Alpha_neg * W)
      
      plot(x = abs(Eta_pos - Eta_neg),
           y = abs(Alpha_pos * W - Alpha_neg * W),
           main = 'Y(-): Eta vs Alpha*W')
      abline(0, 1)
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
      alpha_n = Alpha_neg,
      y_pos_max = max(Y_pos),
      y_neg_max = max(Y_neg),
      y_pos_median = median(Y_pos),
      y_neg_median = median(Y_neg),
      n_tries = counter
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

hist(summary_stats$y_pos_max)
hist(summary_stats$y_neg_max)

summary(summary_stats$y_pos_max)
summary(summary_stats$y_neg_max)

summary(summary_stats$y_pos_median)
summary(summary_stats$y_neg_median)

# Difference in preferential sampling parameters
hist(summary_stats$alpha_p - summary_stats$alpha_n)

# Fractions observed
hist(summary_stats$fraction_observed)
summary(summary_stats$fraction_observed)

write.csv(summary_stats, 'sim_summary_stats.csv')
