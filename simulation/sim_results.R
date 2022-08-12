
library(rstan)
library(aws.s3)
library(stringr)
library(MASS)
library(dplyr)

# ----------------------- Constants -------------------------- #

N_SIMS <- 50

set.seed(123)

FITS_DIR <- 'sim_fits'

DATA_DIR <- 'sim_data'

DO_ALL <- T

KRIGE <- F

# ----------------------- Functions -------------------------- #

# covariance function
cov_exp <- function(dist_matrix, phi, theta){
  
  # Gaussian process covariance matrix
  cov_matrix <- matrix(NA, nrow=nrow(dist_matrix), ncol=ncol(dist_matrix))
  for (i in 1:nrow(dist_matrix)){
    for (j in 1:nrow(dist_matrix)){
      cov_matrix[i,j] <- (phi ** 2) * exp(-dist_matrix[i,j] / theta)
    }
  }
  
  return(cov_matrix)
  
}

krige_gp <- function(gp_samples, theta_samples, phi_samples, d, ids){
  # d: full distance matrix
  # ids: indices of observed data
  # gp samples: gaussian process samples for observed data
  
  # number of cells to extrapolate to
  n_new <- nrow(d) - length(ids)
  
  # number of MCMC samples
  nsamp_old <- nrow(gp_samples)
  
  # array for extrapolated samples
  samples_new <- array(NA, c(nrow(gp_samples), n_new))
  
  cat("generating", nsamp_old, "kriged samples\n")
  progressBar <- txtProgressBar(style = 3)
  percentage.points <- round((1:100/100) * nsamp_old)
  
  for (i in 1:nsamp_old){
    
    theta_i <- theta_samples[i]
    phi_i <- phi_samples[i]
    covmat_i <- cov_exp(d, phi=phi_i, theta=theta_i)
    
    omega_11 <- covmat_i[-ids, -ids]
    omega_12 <- covmat_i[-ids, ids]
    omega_21 <- covmat_i[ids, -ids]
    omega_22 <- covmat_i[ids, ids]
    
    # estimates for observed data
    gp_i <- gp_samples[i,]
    
    # conditional mean and variance
    e_cond <- omega_12 %*% solve(omega_22) %*% gp_i
    var_cond <- omega_11 - omega_12 %*% solve(omega_22) %*% omega_21
    
    # prediction
    gp_pred <- mvrnorm(n=1, mu=e_cond, var_cond)
    samples_new[i,] <- gp_pred
    
    if(i %in% percentage.points){
      setTxtProgressBar(progressBar, i/nsamp_old)
    }
    
  }
  
  return(samples_new)
  
}

merge_krige <- function(gp_obs, gp_krige, ids){
  
  mat_new <- matrix(NA, nrow = nrow(gp_obs), ncol = ncol(gp_obs) + ncol(gp_krige))
  
  mat_new[,ids] <- gp_obs
  
  col_ids <- 1:ncol(mat_new)
  
  ids_predicted <- col_ids[!(col_ids %in% ids)]
  
  mat_new[,ids_predicted] <- gp_krige
  
  return(mat_new)
  
}

# ---------------------- Proposed Model ---------------------- #

fits_slp <- list.files(FITS_DIR)[
  grepl('slp', list.files(FITS_DIR))
  & !(grepl('krige', list.files(FITS_DIR)))
]

print(length(fits_slp))

# Calculate bias and coverage for alpha params
alpha_summary <- c()
for (fit_file in fits_slp){
  
  # get index of fit
  i <- str_split(
    str_split(fit_file, '_')[[1]][3],
    '\\.')[[1]][1]
  
  print(fit_file)
  print(i)
  
  print(paste('simulation ', i))
  
  # Load dataset, params and counts
  params_i <- readRDS(file = paste0(DATA_DIR, '/parameters_', i, '.rds'))
  counts_i <- readRDS(file = paste0(DATA_DIR, '/counts_', i, '.rds'))
  dat_i <- readRDS(file = paste0(DATA_DIR, '/dataset_', i, '.rds'))
  fit_i <- readRDS(paste0(FITS_DIR, '/slp_fit_', i ,'.rds'))
  
  posterior <- as.matrix(fit_i)
  
  # Posterior samples
  alpha_pos_samples <- posterior[,'alpha_pos']
  alpha_neg_samples <- posterior[,'alpha_neg']
  
  # True params
  Alpha_pos <- params_i$Alpha_pos
  Alpha_neg <- params_i$Alpha_neg
  
  bias_pos <- mean(alpha_pos_samples) - Alpha_pos
  bias_neg <- mean(alpha_neg_samples) - Alpha_neg
  
  pos_bounds <- quantile(alpha_pos_samples, c(0.10, 0.90))
  neg_bounds <- quantile(alpha_neg_samples, c(0.10, 0.90))
  
  pos_in_bnds <- 0
  neg_in_bnds <- 0
  if (Alpha_pos > pos_bounds[1] & Alpha_pos < pos_bounds[2]){
    pos_in_bnds <- 1
  }
  if (Alpha_neg > neg_bounds[1] & Alpha_neg < neg_bounds[2]){
    neg_in_bnds <- 1
  }
  
  alpha_summary <- rbind(
    alpha_summary,
    data.frame(
      Bias_pos = bias_pos,
      Bias_neg = bias_neg,
      Pos_covered = pos_in_bnds,
      Neg_covered = neg_in_bnds
    )
  )
  
}
hist(alpha_summary$Bias_pos)
hist(alpha_summary$Bias_neg)
mean(alpha_summary$Pos_covered)
mean(alpha_summary$Neg_covered)
summary(alpha_summary$Bias_pos)
summary(alpha_summary$Bias_neg)

# Now evaluate performance
summary_slp <- c()
for (fit_file in fits_slp){
  
  # get index of fit
  i <- str_split(
    str_split(fit_file, '_')[[1]][3],
    '\\.')[[1]][1]
  
  print(fit_file)
  print(i)
  
  print(paste('simulation ', i))
  
  # Load dataset, params and counts
  params_i <- readRDS(file = paste0(DATA_DIR, '/parameters_', i, '.rds'))
  counts_i <- readRDS(file = paste0(DATA_DIR, '/counts_', i, '.rds'))
  dat_i <- readRDS(file = paste0(DATA_DIR, '/dataset_', i, '.rds'))
  fit_i <- readRDS(paste0(FITS_DIR, '/slp_fit_', i ,'.rds'))
  
  # Extract params
  Theta <- params_i$Theta
  Phi <- params_i$Phi
  Beta_loc <- params_i$Beta_loc
  Beta_pos <- params_i$Beta_pos
  Theta_pos <- params_i$Theta_pos
  Phi_pos <- params_i$Phi_pos
  Beta_neg <- params_i$Beta_neg
  Theta_neg <- params_i$Theta_neg
  Phi_neg <- params_i$Phi_neg
  Alpha_pos <- params_i$Alpha_pos
  Alpha_neg <- params_i$Alpha_neg
  W <- params_i$W
  Eta_pos <- params_i$Eta_pos
  Eta_neg <- params_i$Eta_neg

  posterior <- as.matrix(fit_i)
  
  # Posterior samples
  alpha_pos_samples <- posterior[,'alpha_pos']
  alpha_neg_samples <- posterior[,'alpha_neg']
  w_samples <- posterior[,grepl('w', colnames(posterior))]
  beta_pos_samples <- posterior[,grepl('beta_pos', colnames(posterior))]
  beta_neg_samples <- posterior[,grepl('beta_neg', colnames(posterior))]
  eta_pos_samples <- posterior[,
                               grepl('eta_pos', colnames(posterior))
                               & !(grepl('beta', colnames(posterior)))
                               & !(grepl('theta', colnames(posterior)))]
  eta_neg_samples <- posterior[,
                               grepl('eta_neg', colnames(posterior))
                               & !(grepl('beta', colnames(posterior)))
                               & !(grepl('theta', colnames(posterior)))]
  
  # Posterior means
  alpha_pos_hat <- mean(alpha_pos_samples)
  alpha_neg_hat <- mean(alpha_neg_samples)
  w_hat <- colMeans(w_samples)
  beta_pos_hat <- colMeans(beta_pos_samples)
  beta_neg_hat <- colMeans(beta_neg_samples)
  eta_pos_hat <- colMeans(eta_pos_samples)
  eta_neg_hat <- colMeans(eta_neg_samples)
  
  # check alphas
  par(mfrow=c(1,2))
  hist(posterior[,'alpha_pos'], main=paste('alpha(+)', i))
  abline(v=Alpha_pos, col=2)
  
  hist(posterior[,'alpha_neg'], main=paste('alpha(-)', i))
  abline(v=Alpha_neg, col=2)
  
  summary_i <- data.frame(
    iteration = i,
    alpha_pos = Alpha_pos,
    alpha_neg = Alpha_neg)
  
  if (DO_ALL == T){
    
    if (KRIGE == T){
      
      # check for existing kriged files
      #existing_krige <- list.files(paste0(FITS_DIR, '/krige_slp'))
      #fname <- paste0('eta_pos_all_', i, '.rds')
      existing_krige <- c()
      
      if (!(fname %in% existing_krige)){
      
        # get samples for kriging
        theta_pos_samples <- posterior[,'theta_pos']
        theta_neg_samples <- posterior[,'theta_neg']
        phi_pos_samples <- posterior[,'phi_pos']
        phi_neg_samples <- posterior[,'phi_neg']
      
        # krige disease positive process
        eta_pos_krige <- krige_gp(eta_pos_samples, theta_pos_samples, phi_pos_samples,
                                  dat_i$D, dat_i$I_obs)
        
        # krige disease negative process
        eta_neg_krige <- krige_gp(eta_neg_samples, theta_neg_samples, phi_neg_samples,
                                  dat_i$D, dat_i$I_obs)
      
        eta_pos_samples_all <- merge_krige(eta_pos_samples, eta_pos_krige, dat_i$I_obs)
        eta_neg_samples_all <- merge_krige(eta_neg_samples, eta_neg_krige, dat_i$I_obs)
        
        saveRDS(colMeans(eta_pos_samples_all), paste0(FITS_DIR, '/krige_slp/eta_pos_all_', i, '.rds'))
        saveRDS(colMeans(eta_neg_samples_all), paste0(FITS_DIR, '/krige_slp/eta_neg_all_', i, '.rds'))
      
      } 
      
    }
    
    eta_pos_all <- readRDS(paste0(FITS_DIR, '/krige_slp/eta_pos_all_', i, '.rds'))
    eta_neg_all <- readRDS(paste0(FITS_DIR, '/krige_slp/eta_neg_all_', i, '.rds'))
    
    # calculate log odds
    X <- dat_i$X
    inds <- dat_i$I_obs
    lodds <- (X %*% beta_pos_hat + eta_pos_all + alpha_pos_hat * w_hat) -
      (X %*% beta_neg_hat + eta_neg_all + alpha_neg_hat * w_hat)
    lodds_true <- (X %*% Beta_pos + Eta_pos + Alpha_pos * W) -
      (X %*% Beta_neg + Eta_neg + Alpha_neg * W)
    plot(x=lodds_true, y=lodds, main=paste0('Log Odds, All:', i))
    abline(0, 1, col = 2)
    summary_i$rmse_lodds <- sqrt(mean((lodds - lodds_true) ** 2))
    
    # now calculate log odds at unobserved sites
    lodds_unobs <- lodds[-inds]
    lodds_true_unobs <- lodds_true[-inds]
    summary_i$rmse_lodds_unobs <- sqrt(mean((lodds_unobs - lodds_true_unobs) ** 2))
    
  }
  
  # log odds at observed locations
  X <- dat_i$X
  inds <- dat_i$I_obs
  lodds_obs <- (X[inds,] %*% beta_pos_hat + eta_pos_hat + alpha_pos_hat * w_hat[inds]) -
    (X[inds,] %*% beta_neg_hat + eta_neg_hat + alpha_neg_hat * w_hat[inds])
  lodds_true_obs <- (X[inds,] %*% Beta_pos + Eta_pos[inds] + Alpha_pos * W[inds]) -
    (X[inds,] %*% Beta_neg + Eta_neg[inds] + Alpha_neg * W[inds])
  plot(x=lodds_true_obs, y=lodds_obs, main=paste0('Log Odds, Observed:', i))
  abline(0, 1, col=2)
  summary_i$rmse_lodds_obs <- sqrt(mean((lodds_obs - lodds_true_obs) ** 2))
  
  # compute DIC
  
  # log density of observed data given posterior mean
  mu_pos <- exp(X[inds,] %*% beta_pos_hat + eta_pos_hat + alpha_pos_hat * w_hat[inds])
  mu_neg <- exp(X[inds,] %*% beta_neg_hat + eta_neg_hat + alpha_neg_hat * w_hat[inds])
  log_p_y <- sum(dpois(dat_i$Y_pos, mu_pos, log = TRUE) + dpois(dat_i$Y_neg, mu_neg, log = TRUE))
  
  # compute penalty
  sum_log_p_y <- c()
  for (s in 1:nrow(eta_pos_samples)){

    # vector of means for the ith posterior draw
    lambda_s_pos <- exp(X[inds,] %*% beta_pos_samples[s,] + eta_pos_samples[s,] + alpha_pos_samples[s] * w_samples[s,][inds])
    lambda_s_neg <- exp(X[inds,] %*% beta_neg_samples[s,] + eta_neg_samples[s,] + alpha_neg_samples[s] * w_samples[s,][inds])
    
    # take log density of observed data for the ith sample
    log_p_y_s <- sum(dpois(dat_i$Y_pos, lambda_s_pos, log = T)) + sum(dpois(dat_i$Y_neg, lambda_s_neg, log = T))
    sum_log_p_y <- c(sum_log_p_y, log_p_y_s)
  }
  p_dic <- 2 * (log_p_y - mean(sum_log_p_y))
  
  dic <- -2 * log_p_y + 2 * p_dic
  
  summary_i$dic <- dic
  
  summary_slp <- rbind(
    summary_slp,
    summary_i
  )
  
}

# ---------------------- Poisson Model ----------------------- #

fits_poisson <- list.files(paste0(FITS_DIR, '/', 'poisson'))

summary_poisson <- c()

for (fit_file in fits_poisson){
  
  # get index of fit
  i <- str_split(
    str_split(fit_file, '_')[[1]][3],
    '\\.')[[1]][1]
  
  print(fit_file)
  print(i)
  
  # Load dataset, params and counts
  params_i <- readRDS(file = paste0(DATA_DIR, '/parameters_', i, '.rds'))
  counts_i <- readRDS(file = paste0(DATA_DIR, '/counts_', i, '.rds'))
  dat_i <- readRDS(file = paste0(DATA_DIR, '/dataset_', i, '.rds'))
  fit_i <- readRDS(paste0(FITS_DIR, '/poisson/poisson_fit_', i ,'.rds'))
  
  # Extract params
  Theta <- params_i$Theta
  Phi <- params_i$Phi
  Beta_loc <- params_i$Beta_loc
  Beta_pos <- params_i$Beta_pos
  Theta_pos <- params_i$Theta_pos
  Phi_pos <- params_i$Phi_pos
  Beta_neg <- params_i$Beta_neg
  Theta_neg <- params_i$Theta_neg
  Phi_neg <- params_i$Phi_neg
  Alpha_pos <- params_i$Alpha_pos
  Alpha_neg <- params_i$Alpha_neg
  W <- params_i$W
  Eta_pos <- params_i$Eta_pos
  Eta_neg <- params_i$Eta_neg
  
  # Samples
  posterior <- as.matrix(fit_i)
  beta_pos_samples <- posterior[,grepl('beta_pos', colnames(posterior))]
  beta_neg_samples <- posterior[,grepl('beta_neg', colnames(posterior))]
  
  # Param estimates
  beta_pos_hat <- colMeans(beta_pos_samples)
  beta_neg_hat <- colMeans(beta_neg_samples)
  
  # log odds
  X <- dat_i$X
  lodds <- X %*% beta_pos_hat  - X %*% beta_neg_hat
  lodds_true <- (X %*% Beta_pos + Eta_pos + Alpha_pos * W) -
    (X %*% Beta_neg + Eta_neg + Alpha_neg * W)
  plot(x=lodds_true, y=lodds, main=paste0('Log Odds, All:', i))
  abline(0, 1, col=2)
  
  # at observed locations
  inds <- dat_i$I_obs
  lodds_obs <- lodds[inds]
  lodds_true_obs <- lodds_true[inds]
  plot(x=lodds_true_obs, y=lodds_obs, main=paste0('Log Odds, Observed:', i))
  abline(0, 1, col=2)
  
  # at unobserved locations
  lodds_unobs <- lodds[-inds]
  lodds_true_unobs <- lodds_true[-inds]
  
  # DIC
  
  # log density of observed data given posterior mean
  mu_pos <- exp(X[inds,] %*% beta_pos_hat)
  mu_neg <- exp(X[inds,] %*% beta_neg_hat)
  log_p_y <- sum(dpois(dat_i$Y_pos, mu_pos, log = TRUE) + dpois(dat_i$Y_neg, mu_neg, log = TRUE))
  
  # evaluate penalty
  sum_log_p_y <- c()
  for (s in 1:length(posterior[,'beta_pos[1]'])){
    
    # vector of means for the ith posterior draw
    lambda_s_pos <- exp(X[inds,] %*% beta_pos_samples[s,])
    lambda_s_neg <- exp(X[inds,] %*% beta_neg_samples[s,])
    
    # take log density of observed data for the ith sample
    log_p_y_s <- sum(dpois(dat_i$Y_pos, lambda_s_pos, log = T)) + sum(dpois(dat_i$Y_neg, lambda_s_neg, log = T))
    sum_log_p_y <- c(sum_log_p_y, log_p_y_s)
  }
  p_dic <- 2 * (log_p_y - mean(sum_log_p_y))
  
  dic <- -2 * log_p_y + 2 * p_dic
  
  summary_poisson <- rbind(
    summary_poisson,
    data.frame(
      iteration = i,
      alpha_pos = Alpha_pos,
      alpha_neg = Alpha_neg,
      rmse_lodds = sqrt(mean((lodds - lodds_true) ** 2)),
      rmse_lodds_obs = sqrt(mean((lodds_obs - lodds_true_obs) ** 2)),
      rmse_lodds_unobs = sqrt(mean((lodds_unobs - lodds_true_unobs) ** 2)),
      dic = dic
    )
  )

}

# ---------------------- Spatial Poisson --------------------- #

fits_spois <- list.files(paste0(FITS_DIR, '/', 'spatial_poisson'))

summary_spois <- c()

# Now inspect estimates
for (fit_file in fits_spois){
  
  # get index of fit
  i <- str_split(
    str_split(fit_file, '_')[[1]][4],
    '\\.')[[1]][1]
  
  print(fit_file)
  print(i)
  
  print(paste('simulation ', i))
  
  # Load dataset, params and counts
  params_i <- readRDS(file = paste0(DATA_DIR, '/parameters_', i, '.rds'))
  counts_i <- readRDS(file = paste0(DATA_DIR, '/counts_', i, '.rds'))
  dat_i <- readRDS(file = paste0(DATA_DIR, '/dataset_', i, '.rds'))
  fit_i <- readRDS(paste0(FITS_DIR, '/spatial_poisson/spatial_poisson_fit_', i ,'.rds'))
  
  # Extract params
  Theta <- params_i$Theta
  Phi <- params_i$Phi
  Beta_loc <- params_i$Beta_loc
  Beta_pos <- params_i$Beta_pos
  Theta_pos <- params_i$Theta_pos
  Phi_pos <- params_i$Phi_pos
  Beta_neg <- params_i$Beta_neg
  Theta_neg <- params_i$Theta_neg
  Phi_neg <- params_i$Phi_neg
  Alpha_pos <- params_i$Alpha_pos
  Alpha_neg <- params_i$Alpha_neg
  W <- params_i$W
  Eta_pos <- params_i$Eta_pos
  Eta_neg <- params_i$Eta_neg
  
  # Samples
  posterior <- as.matrix(fit_i)
  beta_pos_samples <- posterior[,grepl('beta_pos', colnames(posterior))]
  beta_neg_samples <- posterior[,grepl('beta_neg', colnames(posterior))]
  eta_pos_samples <- posterior[,
                               grepl('eta_pos', colnames(posterior))
                               & !(grepl('beta', colnames(posterior)))
                               & !(grepl('theta', colnames(posterior)))]
  eta_neg_samples <- posterior[,
                               grepl('eta_neg', colnames(posterior))
                               & !(grepl('beta', colnames(posterior)))
                               & !(grepl('theta', colnames(posterior)))]
  
  # Param estimates
  beta_pos_hat <- colMeans(beta_pos_samples)
  beta_neg_hat <- colMeans(beta_neg_samples)
  eta_pos_hat <- colMeans(eta_pos_samples)
  eta_neg_hat <- colMeans(eta_neg_samples)
  
  summary_i <- data.frame(
    iteration = i,
    alpha_pos = Alpha_pos,
    alpha_neg = Alpha_neg)
  
  if (DO_ALL == T){
    
    if (KRIGE == T){
      
      # check for existing kriged files
      existing_krige <- list.files(paste0(FITS_DIR, '/krige_spatial_poisson'))
      fname <- paste0('eta_pos_all_', i, '.rds')
      
      if (!(fname %in% existing_krige)){
    
        # get samples for kriging
        theta_pos_samples <- posterior[,'theta_pos']
        theta_neg_samples <- posterior[,'theta_neg']
        phi_pos_samples <- posterior[,'phi_pos']
        phi_neg_samples <- posterior[,'phi_neg']
      
        # krige disease positive process
        eta_pos_krige <- krige_gp(eta_pos_samples, theta_pos_samples, phi_pos_samples,
                                  dat_i$D, dat_i$I_obs)
      
        # krige disease negative process
        eta_neg_krige <- krige_gp(eta_neg_samples, theta_neg_samples, phi_neg_samples,
                                  dat_i$D, dat_i$I_obs)
    
        eta_pos_all <- merge_krige(eta_pos_samples, eta_pos_krige, dat_i$I_obs)
        eta_neg_all <- merge_krige(eta_neg_samples, eta_neg_krige, dat_i$I_obs)
        eta_pos_all <- colMeans(eta_pos_all)
        eta_neg_all <- colMeans(eta_neg_all)
        
        saveRDS(eta_pos_all, paste0(FITS_DIR, '/krige_spatial_poisson/eta_pos_all_', i, '.rds'))
        saveRDS(eta_neg_all, paste0(FITS_DIR, '/krige_spatial_poisson/eta_neg_all_', i, '.rds'))
        
      }
    
    } 
      
    eta_pos_all <- readRDS(paste0(FITS_DIR, '/krige_spatial_poisson/eta_pos_all_', i, '.rds'))
    eta_neg_all <- readRDS(paste0(FITS_DIR, '/krige_spatial_poisson/eta_neg_all_', i, '.rds'))
    
    # calculate log odds
    X <- dat_i$X
    inds <- dat_i$I_obs
    lodds <- (X %*% beta_pos_hat + eta_pos_all) -
      (X %*% beta_neg_hat + eta_neg_all)
    lodds_true <- (X %*% Beta_pos + Eta_pos + Alpha_pos * W) -
      (X %*% Beta_neg + Eta_neg + Alpha_neg * W)
    plot(x=lodds_true, y=lodds, main=paste0('Log Odds, All:', i))
    abline(0, 1, col=2)
    summary_i$rmse_lodds <- sqrt(mean((lodds - lodds_true) ** 2))
    
    # now calculate log odds at unobserved sites
    lodds_unobs <- lodds[-inds]
    lodds_true_unobs <- lodds_true[-inds]
    summary_i$rmse_lodds_unobs <- sqrt(mean((lodds_unobs - lodds_true_unobs) ** 2))
    
  }
  
  # log odds at observed locations
  X <- dat_i$X
  inds <- dat_i$I_obs
  lodds_obs <- (X[inds,] %*% beta_pos_hat + eta_pos_hat) -
    (X[inds,] %*% beta_neg_hat + eta_neg_hat)
  lodds_true_obs <- (X[inds,] %*% Beta_pos + Eta_pos[inds] + Alpha_pos * W[inds]) -
    (X[inds,] %*% Beta_neg + Eta_neg[inds] + Alpha_neg * W[inds])
  plot(x=lodds_true_obs, y=lodds_obs, main=paste0('Log Odds, Observed:', i))
  abline(0, 1, col=2)
  summary_i$rmse_lodds_obs = sqrt(mean((lodds_obs - lodds_true_obs) ** 2))
  
  # compute DIC
  
  # log density of observed data given posterior mean
  mu_pos <- exp(X[inds,] %*% beta_pos_hat + eta_pos_hat)
  mu_neg <- exp(X[inds,] %*% beta_neg_hat + eta_neg_hat)
  log_p_y <- sum(dpois(dat_i$Y_pos, mu_pos, log = TRUE) + dpois(dat_i$Y_neg, mu_neg, log = TRUE))
  
  # evaluate penalty
  sum_log_p_y <- c()
  for (s in 1:nrow(eta_pos_samples)){
    
    # vector of means for the ith posterior draw
    lambda_s_pos <- exp(X[inds,] %*% beta_pos_samples[s,] + eta_pos_samples[s,])
    lambda_s_neg <- exp(X[inds,] %*% beta_neg_samples[s,] + eta_neg_samples[s,])
    
    # take log density of observed data for the ith sample
    log_p_y_s <- sum(dpois(dat_i$Y_pos, lambda_s_pos, log = T)) + sum(dpois(dat_i$Y_neg, lambda_s_neg, log = T))
    sum_log_p_y <- c(sum_log_p_y, log_p_y_s)
  }
  p_dic <- 2 * (log_p_y - mean(sum_log_p_y))
  
  dic <- -2 * log_p_y + 2 * p_dic
  
  summary_i$dic <- dic
  
  summary_spois <- rbind(
    summary_spois,
    summary_i
  )
  
}

# ----------------- DIC Comparison --------------------------- #

summary_spois$model <- 'Spatial Poisson'
summary_poisson$model <- 'Poisson'
summary_slp$model <- 'SLP'

dat_combined <- rbind(summary_spois,
                      summary_poisson,
                      summary_slp)

ggplot(dat_combined, aes(x=model, y=dic)) + 
  geom_boxplot()

dat_combined %>%
  group_by(model) %>%
  summarise(mean_dic = mean(dic),
            median_dic = median(dic),
            q25_dic = quantile(dic, 0.25),
            q75_dic = quantile(dic, 0.75))
# model           mean_dic median_dic q25_dic q75_dic
 # Poisson           30593.     15892.   5529.  31150.
 # SLP                 929.       939.    799.   1106.
 # Spatial Poisson     957.       968.    835.   1116.

# fraction of times each model has lowest DIC
slp_lowest <- c()
poisson_lowest <- c()
spois_lowest <- c()
for (i in unique(dat_combined$iteration)){
  
  slp_dic <- dat_combined %>%
    filter(model == 'SLP',
           iteration == i) %>%
    pull(dic)
  poisson_dic <- dat_combined %>%
    filter(model == 'Poisson',
           iteration == i) %>%
    pull(dic)
  spatial_poisson_dic <- dat_combined %>%
    filter(model == 'Spatial Poisson',
           iteration == i) %>%
    pull(dic)
  
  if (min(c(slp_dic, poisson_dic, spatial_poisson_dic)) == slp_dic){
    slp_lowest <- c(slp_lowest, 1)
    poisson_lowest <- c(poisson_lowest, 0)
    spois_lowest <- c(spois_lowest, 0)
  } else if (min(c(slp_dic, poisson_dic, spatial_poisson_dic)) == spatial_poisson_dic){
    slp_lowest <- c(slp_lowest, 0)
    poisson_lowest <- c(poisson_lowest, 0)
    spois_lowest <- c(spois_lowest, 1)
  }
  
  
  
}

# pairwise differences in DIC for each iteration
slp_spois_diffs <- c()
slp_pois_diffs <- c()
spois_pois_diffs <- c()

for (i in unique(dat_combined$iteration)){
  
  slp_dic <- dat_combined %>%
    filter(model == 'SLP',
           iteration == i) %>%
    pull(dic)
  poisson_dic <- dat_combined %>%
    filter(model == 'Poisson',
           iteration == i) %>%
    pull(dic)
  spatial_poisson_dic <- dat_combined %>%
    filter(model == 'Spatial Poisson',
           iteration == i) %>%
    pull(dic)
  
  slp_spois_diffs <- c(
    slp_spois_diffs, slp_dic - spatial_poisson_dic
  )
  slp_pois_diffs <- c(
    slp_pois_diffs, slp_dic - poisson_dic
  )
  spois_pois_diffs <- c(
    spois_pois_diffs, spatial_poisson_dic - poisson_dic
  )
  
}

print(summary(slp_spois_diffs))
print(summary(slp_spois_diffs))

# ------------------------ Compare RMSE ------------------------ #

# Overall RMSE
dat_combined %>%
  group_by(model) %>%
  summarise(mean_rmse = mean(rmse_lodds),
            median_rmse = median(rmse_lodds),
            sd = sd(rmse_lodds))

dat_combined <- dat_combined %>%
  mutate(abs_alpha_diff = abs(alpha_pos - alpha_neg))

# RMSE diff low
dat_combined %>%
  filter(abs_alpha_diff < quantile(abs_alpha_diff, 0.33)) %>%
  group_by(model) %>%
  summarise(mean_rmse_un = mean(rmse_lodds_unobs),
            median_rmse_un = median(rmse_lodds_unobs),
            sd_un = sd(rmse_lodds_unobs))

# RMSE diff med
dat_combined %>%
  filter(abs_alpha_diff > quantile(abs_alpha_diff, 0.33),
         abs_alpha_diff < quantile(abs_alpha_diff, 0.66)) %>%
  group_by(model) %>%
  summarise(mean_rmse_un = mean(rmse_lodds_unobs),
            median_rmse_un = median(rmse_lodds_unobs),
            sd_un = sd(rmse_lodds_unobs))

# RMSE diff high
dat_combined %>%
  filter(abs_alpha_diff > quantile(abs_alpha_diff, 0.66)) %>%
  group_by(model) %>%
  summarise(mean_rmse_un = mean(rmse_lodds_unobs),
            median_rmse_un = median(rmse_lodds_unobs),
            sd_un = sd(rmse_lodds_unobs))

dat_combined <- dat_combined %>%
  mutate(model = case_when(
    model == "Poisson" ~ 'Pois-R',
    model == "Spatial Poisson" ~ 'Spat-P',
    TRUE ~ 'Pref-S'
  ))

dat_combined$Model <- factor(dat_combined$model, levels = c('Pois-R', 
                                                            'Spat-P',
                                                            'Pref-S'))

# Lineplot
ggplot(data=dat_combined, 
       aes(x=abs_alpha_diff, 
           y=rmse_lodds, 
           group=Model, 
           colour=Model)) +
  geom_line() +
  geom_point() +
  labs(x = 'Absolute Difference',
       y = "RMSE") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

dat_all <- dat_combined %>%
  mutate(Data = 'All')
dat_all <- rbind(dat_all,
                 dat_combined %>%
                 mutate(Data = case_when(
                   abs_alpha_diff < quantile(abs_alpha_diff, 0.33) ~ 'Low',
                   abs_alpha_diff > quantile(abs_alpha_diff, 0.66) ~ 'High',
                   TRUE ~ 'Med'
                 ))
                 )


# Boxplot
ggplot() + 
  geom_boxplot(
    dat_all %>% filter(Data != 'Med'), 
    mapping = aes(Data, rmse_lodds, fill=Model)) + 
  scale_x_discrete(limits=c("All", "Low", "High")) + 
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y='RMSE')

dat_combined %>%
  mutate()

# Table summarizing metrics
summary_tab = dat_all %>% 
  filter(Sampling != 'Med') %>%
  group_by(Sampling, Model) %>%
  summarise(
    mean_rmse = round(mean(rmse_lodds), 3),
    sd_rmse = round(sd(rmse_lodds), 3),
    mean_dic = round(mean(dic), 3),
    sd_dic = round(sd(dic), 3)
  )
