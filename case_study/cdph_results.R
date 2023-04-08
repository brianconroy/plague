
# This script generates risk maps from the models
# fit to the CDPH plague surveillance data.

library(rstan)
library(aws.s3)
library(sp)
library(raster)
library(MASS)
library(parallel)
library(doParallel)
library(tictoc)
library(fields)
library(np)
library(ggplot2)
library(bayesplot)
library(posterior)

# ----------------------- Constants ------------------------- #

set.seed(123)

AGG_FACTOR <- 7

DOWNLOAD_FROM_S3 <- TRUE

SLP_TAG <- 'slp' 
SPOIS_TAG <- 'spois'
POISSON_TAG <- 'pois'

DATA_FILE <- 'data/plague_data_agg7.rds'

# STAN fits
POISSON_FILE <- paste0('cdph_fits/plague_poisson_fit_', POISSON_TAG, '.rds')
SLP_FILE <- paste0('cdph_fits/plague_slp_fit_', SLP_TAG, '.rds')
SPATIAL_POISSON_FILE <- paste0('cdph_fits/plague_spatial_poisson_fit_', SPOIS_TAG, '.rds')

ANALYZE_POIS <- TRUE
ANALYZE_SPOIS <- TRUE
ANALYZE_SLP <- TRUE

# Kriging
SLP_KRIGE_POS_NAME <- paste0('cdph_fits/cdph_slp_krige_pos_', SLP_TAG, '.rds')
SLP_KRIGE_NEG_NAME <- paste0('cdph_fits/cdph_slp_krige_neg_', SLP_TAG, '.rds')
SPATIAL_POISSON_KRIGE_POS_NAME <- paste0('cdph_fits/cdph_spois_krige_pos_', SPOIS_TAG, '.rds')
SPATIAL_POISSON_KRIGE_NEG_NAME <- paste0('cdph_fits/cdph_spois_krige_neg_', SPOIS_TAG, '.rds')

KRIGE_SPOIS <- TRUE
KRIGE_SLP <- TRUE

INTERPOLATE_SPOIS <- TRUE
INTERPOLATE_SLP <- TRUE

FIND_BWS_SPOIS <- TRUE
FIND_BWS_SLP <- TRUE

# Outputs
SLP_FINE_RASTER_FILE <- paste0('cdph_fits/cdph_slp_raster_fine_', SLP_TAG, '.rds')
SLP_SIGNIFICANCE_RASTER_FILE <- paste0('cdph_fits/cdph_slp_raster_significance_', SLP_TAG, '.rds')

SPOIS_FINE_RASTER_FILE <- paste0('cdph_fits/cdph_spois_raster_', SPOIS_TAG, '.rds')
SPOIS_SIGNIFICANCE_RASTER_FILE <- paste0('cdph_fits/cdph_spois_raster_significance_', SPOIS_TAG, '.rds')

POISSON_FINE_RASTER_FILE <- 'cdph_fits/cdph_poisson_raster_fine.rds'
POISSON_SIGNIFICANCE_RASTER_FILE <- 'cdph_fits/cdph_poisson_raster_significance.rds'

if (DOWNLOAD_FROM_S3){
  source('keys.R')
  Sys.setenv(
    AWS_ACCESS_KEY_ID = AWS_ACCESS_KEY_ID,
    AWS_SECRET_ACCESS_KEY = AWS_SECRET_ACCESS_KEY,
    AWS_REGION = AWS_REGION
  )
}

# ----------------------- Functions ------------------------ #

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
    gp_pred <- mvrnorm(n=1, mu=e_cond, Sigma=var_cond)
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

get_optimal_bws <- function(r_pred, r_train, z){
  
  txdat <- data.frame(xyFromCell(r_train, (1:ncell(r_train))[!is.na(r_train[])]))
  x <- txdat[,1]
  y <- txdat[,2]
  df_new <- data.frame(xyFromCell(r_pred, (1:ncell(r_pred))[!is.na(r_pred[])]))
  
  bws_mod <- npregbw(
    formula = z~x+y,
    regtype = "lc",
    ckertype = "gaussian")
  
  return(bws_mod$bw)
  
}

interpolate_gp_batched  <- function(samples, bws, r_train, r_pred, batch_size = 500){
  
  samples_int <- lapply(1:(nrow(samples)/batch_size),function(batch) {
    print(paste(batch))
    interpolate_gp(samples[(batch_size * (batch-1) + 1):(batch_size * batch),], bws, r_train, r_pred)
  })
  samples_int <- do.call(rbind, samples_int)
  
  return(samples_int)
  
}

interpolate_gp <- function(samples, bws, r_train, r_pred, out_file = NULL){
  
  txdat <- data.frame(xyFromCell(r_train, (1:ncell(r_train))[!is.na(r_train[])]))
  x <- txdat[,1]
  y <- txdat[,2]
  df_new <- data.frame(xyFromCell(r_pred, (1:ncell(r_pred))[!is.na(r_pred[])]))
  
  cores <- detectCores()
  cl <- makeCluster(cores[1]-1)
  registerDoParallel(cl)
  
  print(Sys.time())
  print(paste("interpolating", nrow(samples), "samples"))
  tic("smoothing")
  samples_pred <- foreach(i=1:nrow(samples), .combine=rbind) %dopar% {
    library(np)
    z <- samples[i,]
    model.np <- npreg(bws=bws,
                      formula=z~x+y,
                      regtype="lc",
                      ckertype="gaussian")
    pred <- predict(model.np,
                    newdata=df_new)
    pred
  }
  toc()
  
  stopCluster(cl)
  
  if (!is.null(out_file)){
    save_output(samples_pred, out_file)
  }
  
  return(samples_pred)
  
}

calc_posterior_risk <- function(x, beta_pos_samp, beta_neg_samp, eta_pos_samp,
                                eta_neg_samp, alpha_pos_samp, alpha_neg_samp,
                                w_samp){
  
  n.samp <- nrow(beta_pos_samp)
  n.cell <- ncol(w_samp)
  risk_samp <- array(NA, c(n.samp, n.cell))
  
  for (i in 1:n.samp){
    
    beta_ca <- beta_pos_samp[i,]
    beta_co <- beta_neg_samp[i,]
    alpha_ca <- alpha_pos_samp[i]
    alpha_co <- alpha_neg_samp[i]
    eta_ca <- eta_pos_samp[i,]
    eta_co <- eta_neg_samp[i,]
    w <- w_samp[i,]
    
    lodds.i <- (x %*% beta_ca + eta_ca + alpha_ca * w) - 
      (x %*% beta_co + eta_co + alpha_co * w)
    risk_samp[i,] <- t(exp(lodds.i)/(1 + exp(lodds.i)))
    
  }
  
  return(risk_samp)
  
}

write_latex_table <- function(df, fname, path){
  
  df[,] <- lapply(df[, ], as.character)
  tab <- "\\begin{table}[h]
  \\begin{center}
  \\begin{tabular}{l*{5}{c}r}\n"
  header <- paste(names(df), collapse=" & ")
  tab <- paste(tab, header, "\\\\ \n \\hline \n")
  for (i in 1:nrow(df)){
    row_ <- paste(df[i,], collapse=" & ")
    row_ <- paste(row_, '\\\\\n')
    tab <- paste(tab, row_, sep=" ")
  }
  end <- "\\hline
  \\end{tabular}
  \\caption[placeholder]
  {
  placeholder
  }
  \\label{tab:placeholder}
  \\end{center}
  \\end{table}"
  tab <- paste(tab, end, sep=" ")
  path <- paste(path, fname, sep="")
  write(tab, path)
  
}

calc_sig_map <- function(risk_df, threshold){
  
  r_inds <- apply(risk_df, 2, function(x) return(mean(x > threshold) > 0.95))
  r_inds_raster <- caPr[[1]]
  r_inds_raster[][!is.na(r_inds_raster[])] <- r_inds
  return(r_inds_raster)
  
}

# --------------------- Covariates ------------------------- #

# Spatial covariates from PRISM
caPr <- raster::stack("data/prism_pcas_ca.grd")
caPr.disc <- aggregate(caPr, fact = AGG_FACTOR)
plot(caPr.disc)

N <- length(caPr.disc[[1]][][!is.na(caPr.disc[[1]][])])
print(N)

##############################################################
#                       Create Maps                          #
##############################################################

model_comparison <- c()

# --------------------- Poisson Model ---------------------- #

if (ANALYZE_POIS){
  
  dat <- readRDS(file = DATA_FILE)
  fit_pois <- readRDS(file = POISSON_FILE)
  
  print(paste('N:', dat$N))
  
  posterior <- as.matrix(fit_pois)
  beta_pos_hat <- colMeans(posterior[,grepl('beta_pos', colnames(posterior))])
  beta_neg_hat <- colMeans(posterior[,grepl('beta_neg', colnames(posterior))])
  
  # Posterior risk samples
  X_f <- cbind(1, 
               caPr[[1]][][!is.na(caPr[[1]][])],
               caPr[[2]][][!is.na(caPr[[2]][])])
  beta_pos_samp <- posterior[,grepl('beta_pos', colnames(posterior))]
  beta_neg_samp <- posterior[,grepl('beta_neg', colnames(posterior))]
  pois_risk_samp <- calc_posterior_risk(X_f, 
                                        beta_pos_samp, 
                                        beta_neg_samp, 
                                        eta_pos_samp = matrix(0, ncol = nrow(X_f), nrow = nrow(beta_pos_samp)),
                                        eta_neg_samp = matrix(0, ncol = nrow(X_f), nrow = nrow(beta_pos_samp)),
                                        alpha_pos_samp = matrix(0, nrow = nrow(X_f) , ncol = 1), 
                                        alpha_neg_samp = matrix(0, nrow = nrow(X_f) , ncol = 1),
                                        w_samp = matrix(0, nrow = nrow(beta_pos_samp), ncol = nrow(X_f)))
  
  # Risk map
  r_pois_f <- caPr[[1]]
  r_pois_f[][!is.na(r_pois_f[])] <- colMeans(pois_risk_samp)
  plot(r_pois_f)

  # Significance indicator maps  
  r_pois_ind1 <- calc_sig_map(pois_risk_samp, threshold = 0.064)
  r_pois_ind2 <- calc_sig_map(pois_risk_samp, threshold = 0.032)
  r_pois_ind <- stack(r_pois_ind1, r_pois_ind2)
  plot(r_pois_ind)

  saveRDS(r_pois_f, POISSON_FINE_RASTER_FILE)
  saveRDS(r_pois_ind, POISSON_SIGNIFICANCE_RASTER_FILE)

  # DIC
  
  X <- dat$X
  
  # Offset
  N_offset <- dat$Y_loc[dat$I_obs]
  
  # Log density of observed data given posterior mean
  mu_pos <- N_offset * exp(X[dat$I_obs,] %*% beta_pos_hat)
  mu_neg <- N_offset * exp(X[dat$I_obs,] %*% beta_neg_hat)
  log_p_y <- sum(dpois(dat$Y_pos, mu_pos, log = TRUE) + dpois(dat$Y_neg, mu_neg, log = TRUE))
  
  # Evaluate penalty
  sum_log_p_y <- c()
  for (s in 1:length(posterior[,'beta_pos[1]'])){
    
    # Vector of means for the ith posterior draw
    lambda_s_pos <- N_offset * exp(X[dat$I_obs,] %*% beta_pos_samp[s,])
    lambda_s_neg <- N_offset * exp(X[dat$I_obs,] %*% beta_neg_samp[s,])
    
    # Take log density of observed data for the ith sample
    log_p_y_s <- sum(dpois(dat$Y_pos, lambda_s_pos, log = T)) + sum(dpois(dat$Y_neg, lambda_s_neg, log = T))
    sum_log_p_y <- c(sum_log_p_y, log_p_y_s)
    
  }
  p_dic <- 2 * (log_p_y - mean(sum_log_p_y))
  
  dic <- -2 * log_p_y + 2 * p_dic
  
  model_comparison <- rbind(
    model_comparison,
    data.frame(
      model = 'Poisson',
      dic = dic
    )
  )
  
  # Parameter estimates
  pois_params <- c()
  for (i in 1:3){
    pois_params <- rbind(
      pois_params,
      data.frame(
        Parameter = paste0('Beta (+,', i, ')'),
        Posterior.Mean = round(unname(beta_pos_hat[i]), 3),
        Posterior.80.lb = round(unname(quantile(beta_pos_samp[,i], 0.10)), 3),
        Posterior.80.ub = round(unname(quantile(beta_pos_samp[,i], 0.90), 3))
      )
    )
  }
  for (i in 1:3){
    pois_params <- rbind(
      pois_params,
      data.frame(
        Parameter = paste0('Beta (-,', i, ')'),
        Posterior.Mean = round(unname(beta_neg_hat[i]), 3),
        Posterior.80.lb = round(unname(quantile(beta_neg_samp[,i], 0.10)), 3),
        Posterior.80.ub = round(unname(quantile(beta_neg_samp[,i], 0.90)), 3)
      )
    )
  }
  write_latex_table(pois_params, '/poisson_estimates.txt', getwd())
  
  # Other checks
  X <- dat$X
  par(mfrow = c(1, 2))
  exp_pos <- exp(X[dat$I_obs,] %*% beta_pos_hat)
  plot(x = dat$Y_pos, y = exp_pos, xlab = 'Observed Count', ylab = 'Expected Count',
       main = 'Poisson (Disease Positive)')
  abline(0, 1, col = 2)
  
  exp_neg <- exp(X[dat$I_obs,] %*% beta_neg_hat)
  plot(x = dat$Y_neg, y = exp_neg, xlab = 'Observed Count', ylab = 'Expected Count',
       main = 'Poisson (Disease Negative)')
  abline(0, 1, col = 2)
  
  plot_title <- ggtitle("Posterior distributions",
                        "with medians and 80% intervals")
  color_scheme_set("green")
  mcmc_areas(posterior,
             pars = c("beta_pos[1]", "beta_pos[2]", "beta_pos[3]"),
             prob = 0.8) + plot_title
  
  plot_title <- ggtitle("Posterior distributions",
                        "with medians and 80% intervals")
  color_scheme_set("green")
  mcmc_areas(posterior,
             pars = c("beta_neg[1]", "beta_neg[2]", "beta_neg[3]"),
             prob = 0.8) + plot_title
  
}

# ----------------- Spatial Poisson Model ------------------ #

if (ANALYZE_SPOIS){

  if (DOWNLOAD_FROM_S3){
  
    # Retrieve fit
    save_object(
      bucket = "plague-analysis",
      object = gsub('cdph_fits/', '', SPATIAL_POISSON_FILE),
      file = SPATIAL_POISSON_FILE
    )

  }
  
  dat <- readRDS(file = DATA_FILE)
  fit_spois <- readRDS(file = SPATIAL_POISSON_FILE)
  
  print(paste('N:', dat$N))
  print(paste('N_obs:', dat$N_obs))
  
  # Get estimates
  posterior <- as.matrix(fit_spois)
  beta_pos_samp <- posterior[,grepl('beta_pos', colnames(posterior))]
  beta_neg_samp <- posterior[,grepl('beta_neg', colnames(posterior))]
  beta_pos_hat <- colMeans(beta_pos_samp)
  beta_neg_hat <- colMeans(beta_neg_samp)
  eta_pos_samples <- posterior[,
                               grepl('eta_pos', colnames(posterior))
                               & !(grepl('beta', colnames(posterior)))
                               & !(grepl('theta', colnames(posterior)))]
  eta_pos_hat <- colMeans(eta_pos_samples)
  eta_neg_samples <- posterior[,
                               grepl('eta_neg', colnames(posterior))
                               & !(grepl('beta', colnames(posterior)))
                               & !(grepl('theta', colnames(posterior)))]
  eta_neg_hat <- colMeans(eta_neg_samples)
  beta_pos_samples <- posterior[,grepl('beta_pos', colnames(posterior))]
  beta_neg_samples <- posterior[,grepl('beta_neg', colnames(posterior))]
  
  print(paste('MCMC draws:',length(posterior[,'theta_pos'])))
  
  hist(eta_pos_hat)
  hist(eta_neg_hat)
  
  # DIC
  
  # Offset
  N_offset <- dat$Y_loc[dat$I_obs]
  
  # Log density of observed data given posterior mean
  X <- dat$X
  inds <- dat$I_obs
  mu_pos <- N_offset * exp(X[inds,] %*% beta_pos_hat + eta_pos_hat)
  mu_neg <- N_offset * exp(X[inds,] %*% beta_neg_hat + eta_neg_hat)
  log_p_y <- sum(dpois(dat$Y_pos, mu_pos, log = TRUE) + dpois(dat$Y_neg, mu_neg, log = TRUE))
  
  # Evaluate penalty
  sum_log_p_y <- c()
  for (s in 1:nrow(eta_pos_samples)){
    
    # Vector of means for the ith posterior draw
    lambda_s_pos <- N_offset * exp(X[inds,] %*% beta_pos_samp[s,] + eta_pos_samples[s,])
    lambda_s_neg <- N_offset * exp(X[inds,] %*% beta_neg_samp[s,] + eta_neg_samples[s,])
    
    # Take log density of observed data for the ith sample
    log_p_y_s <- sum(dpois(dat$Y_pos, lambda_s_pos, log = T)) + sum(dpois(dat$Y_neg, lambda_s_neg, log = T))
    sum_log_p_y <- c(sum_log_p_y, log_p_y_s)
  }
  p_dic <- 2 * (log_p_y - mean(sum_log_p_y))
  
  dic <- -2 * log_p_y + 2 * p_dic
  
  model_comparison <- rbind(
    model_comparison,
    data.frame(
      model = 'Spatial Poisson',
      dic = dic
    )
  )
  
  if (KRIGE_SPOIS){
  
    # Krige eta random effects
    theta_pos_samples <- posterior[,'theta_pos']
    theta_neg_samples <- posterior[,'theta_neg']
    phi_pos_samples <- posterior[,'phi_pos']
    phi_neg_samples <- posterior[,'phi_neg']
    
    hist(theta_pos_samples)
    hist(theta_neg_samples)
    hist(phi_pos_samples)
    hist(phi_neg_samples)
    
    eta_neg_krige <- krige_gp(eta_neg_samples,
                              theta_samples = theta_neg_samples,
                              phi_samples = phi_neg_samples,
                              d = dat$D,
                              ids = dat$I_obs)
    
    eta_pos_krige <- krige_gp(eta_pos_samples,
                              theta_samples = theta_pos_samples,
                              phi_samples = phi_pos_samples,
                              d = dat$D,
                              ids = dat$I_obs)
  
  }
  
  # Combine observed and kriged samples
  eta_neg_all <- merge_krige(eta_neg_samples, eta_neg_krige, ids = dat$I_obs)
  eta_pos_all <- merge_krige(eta_pos_samples, eta_pos_krige, ids = dat$I_obs)
  
  if (INTERPOLATE_SPOIS){
    
    # Define a raster at the desired resolution
    r_pred <- caPr[[1]]
    
    # Define a raster at the starting (coarse) res.
    r_train <- aggregate(r_pred, fact = AGG_FACTOR)
    
    # Get optimal bandwidthts for eta (+), eta (-), and w
    bws_eta_neg <- get_optimal_bws(r_pred, r_train, colMeans(eta_neg_all))
    bws_eta_pos <- get_optimal_bws(r_pred, r_train, colMeans(eta_pos_all))

    # Interpolate eta (+) samples
    eta_pos_interp <- interpolate_gp_batched(eta_pos_all, bws_eta_pos, r_train, r_pred, batch_size = 500)
    
    # Interpolate eta (-) samples
    eta_neg_interp <- interpolate_gp_batched(eta_neg_all, bws_eta_neg, r_train, r_pred, batch_size = 500)
    
  }

  # Plague risk map - fine resolution
  X_f <- cbind(1,
               caPr[[1]][][!is.na(caPr[[1]][])],
               caPr[[2]][][!is.na(caPr[[2]][])])
  risk_spois_f <- calc_posterior_risk(X_f, beta_pos_samples, beta_neg_samples, eta_pos_interp,
                                      eta_neg_interp,
                                      alpha_pos_samp = matrix(0, nrow = nrow(X_f) , ncol = 1), 
                                      alpha_neg_samp = matrix(0, nrow = nrow(X_f) , ncol = 1),
                                      w_samp = matrix(0, nrow = nrow(X_f), ncol = ncol(eta_neg_interp)))
  
  r_spois_f <- caPr[[1]]
  r_spois_f[][!is.na(r_spois_f[])] <- colMeans(risk_spois_f)
  plot(r_spois_f)
  
  # Significance indicators
  r_spois_ind1 <- calc_sig_map(risk_spois_f, threshold = 0.064)
  r_spois_ind2 <- calc_sig_map(risk_spois_f, threshold = 0.032)
  r_spois_ind <- stack(r_spois_ind1, r_spois_ind2)
  plot(r_spois_ind)
  
  saveRDS(r_spois_f, SPOIS_FINE_RASTER_FILE)
  saveRDS(r_spois_ind, SPOIS_SIGNIFICANCE_RASTER_FILE)

  spois_params <- c()
  for (i in 1:3){
    spois_params <- rbind(
      spois_params,
      data.frame(
        Parameter = paste0('Beta (+,', i, ')'),
        Posterior.Mean = round(unname(beta_pos_hat[i]), 3),
        Posterior.80.lb = round(unname(quantile(beta_pos_samples[,i], 0.10)), 3),
        Posterior.80.ub = round(unname(quantile(beta_pos_samples[,i], 0.90), 3))
      )
    )
  }
  for (i in 1:3){
    spois_params <- rbind(
      spois_params,
      data.frame(
        Parameter = paste0('Beta (-,', i, ')'),
        Posterior.Mean = round(unname(beta_neg_hat[i]), 3),
        Posterior.80.lb = round(unname(quantile(beta_neg_samples[,i], 0.10)), 3),
        Posterior.80.ub = round(unname(quantile(beta_neg_samples[,i], 0.90)), 3)
      )
    )
  }
  write_latex_table(spois_params, '/spois_estimates.txt', getwd())
  
  # Other checks
  plot_title <- ggtitle("Posterior distributions",
                        "with medians and 80% intervals")
  color_scheme_set("green")
  mcmc_areas(posterior,
             pars = c("alpha_pos", "alpha_neg"),
             prob = 0.8) + plot_title
  
}

# --------------------- Proposed Model --------------------- #

if (ANALYZE_SLP){
  
  if (DOWNLOAD_FROM_S3){
    
    # retrieve fit
    save_object(
      bucket = "plague-analysis",
      object = gsub('cdph_fits/', '', SLP_FILE),
      file = SLP_FILE
    )
    
  }
  
  dat <- readRDS(file = DATA_FILE)
  fit_slp <- readRDS(file = SLP_FILE)
  
  print(dat$N)
  
  # Get posterior samples
  posterior <- as.matrix(fit_slp)
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
  w_hat <- colMeans(w_samples)
  
  hist(alpha_pos_samples)
  hist(alpha_neg_samples)
  
  # DIC
  
  # Offset
  N_offset <- dat$Y_loc[dat$I_obs]
  
  X <- dat$X
  inds <- dat$I_obs
  # Log density of observed data given posterior mean
  mu_pos <- N_offset * exp(X[inds,] %*% beta_pos_hat + eta_pos_hat + alpha_pos_hat * w_hat[inds])
  mu_neg <- N_offset * exp(X[inds,] %*% beta_neg_hat + eta_neg_hat + alpha_neg_hat * w_hat[inds])
  log_p_y <- sum(dpois(dat$Y_pos, mu_pos, log = TRUE) + dpois(dat$Y_neg, mu_neg, log = TRUE))
  
  # Compute penalty
  sum_log_p_y <- c()
  for (s in 1:nrow(eta_pos_samples)){
    
    # Vector of means for the ith posterior draw
    lambda_s_pos <- N_offset * exp(X[inds,] %*% beta_pos_samples[s,] + eta_pos_samples[s,] + alpha_pos_samples[s] * w_samples[s,][inds])
    lambda_s_neg <- N_offset * exp(X[inds,] %*% beta_neg_samples[s,] + eta_neg_samples[s,] + alpha_neg_samples[s] * w_samples[s,][inds])
    
    # Take log density of observed data for the ith sample
    log_p_y_s <- sum(dpois(dat$Y_pos, lambda_s_pos, log = T)) + sum(dpois(dat$Y_neg, lambda_s_neg, log = T))
    sum_log_p_y <- c(sum_log_p_y, log_p_y_s)
  }
  p_dic <- 2 * (log_p_y - mean(sum_log_p_y))
  
  dic <- -2 * log_p_y + 2 * p_dic
  
  # Credible interval for alpha(-) - alpha(-)
  alpha_diff <- alpha_pos_samples - alpha_neg_samples
  alpha_diff_est <- list(
    mu = round(mean(alpha_diff), 3),
    ci_bounds = round(quantile(alpha_diff, c(0.025, 0.975)), 3)
  )
  
  model_comparison <- rbind(
    model_comparison,
    data.frame(
      model = 'SLP',
      dic = dic
    )
  )
  
  hist(w_samples)
  hist(eta_pos_samples)
  hist(eta_neg_samples)
 
  if (KRIGE_SLP){
    
    # Krige eta random effects
    theta_pos_samples <- posterior[,'theta_pos']
    theta_neg_samples <- posterior[,'theta_neg']
    phi_pos_samples <- posterior[,'phi_pos']
    phi_neg_samples <- posterior[,'phi_neg']
    
    hist(theta_pos_samples)
    hist(theta_neg_samples)
    hist(phi_pos_samples)
    hist(phi_neg_samples)
    
    hist(posterior[,'theta'])
    hist(posterior[,'phi'])
    
    eta_neg_krige <- krige_gp(eta_neg_samples,
                              theta_samples = theta_neg_samples,
                              phi_samples = phi_neg_samples,
                              d = dat$D,
                              ids = dat$I_obs)
    
    eta_pos_krige <- krige_gp(eta_pos_samples,
                              theta_samples = theta_pos_samples,
                              phi_samples = phi_pos_samples,
                              d = dat$D,
                              ids = dat$I_obs)
    
  }
  
  # Combine observed and kriged samples
  eta_neg_all <- merge_krige(eta_neg_samples, eta_neg_krige, ids = dat$I_obs)
  eta_pos_all <- merge_krige(eta_pos_samples, eta_pos_krige, ids = dat$I_obs)
  
  hist(eta_neg_all)
  hist(eta_pos_all)
  hist(colMeans(w_samples))
  
  if (INTERPOLATE_SLP){
    
    # Define a raster at the desired resolution
    r_pred <- caPr[[1]]
    
    # Define a raster at the starting (coarse) res.
    r_train <- aggregate(r_pred, fact=AGG_FACTOR)
    
    # Get optimal bandwidthts for eta (+), eta (-), and w
    bws_eta_neg <- get_optimal_bws(r_pred, r_train, colMeans(eta_neg_all))
    bws_eta_pos <- get_optimal_bws(r_pred, r_train, colMeans(eta_pos_all))
    bws_w <- get_optimal_bws(r_pred, r_train, colMeans(w_samples))
    
    # Interpolate each posterior sample of w via 2d kernel smoothing
    w_interp <- interpolate_gp_batched(w_samples, bws_w, r_train, r_pred, batch_size = 500)

    # Interpolate eta (+)
    eta_pos_interp <- interpolate_gp_batched(eta_pos_all, bws_eta_pos, r_train, r_pred, batch_size = 500)

    # Interpolate eta (-)
    eta_neg_interp <- interpolate_gp_batched(eta_neg_all, bws_eta_neg, r_train, r_pred, batch_size = 500)

  }
  
  # Calculate risk
  X_f <- cbind(1,
               caPr[[1]][][!is.na(caPr[[1]][])],
               caPr[[2]][][!is.na(caPr[[2]][])])
  risk_slp_f <- calc_posterior_risk(X_f, beta_pos_samples, beta_neg_samples, eta_pos_interp,
                                    eta_neg_interp, alpha_pos_samples, alpha_neg_samples, w_interp)
  
  # Risk map
  r_slp_f <-  caPr[[1]]
  r_slp_f[!is.na(r_slp_f[])] <- colMeans(risk_slp_f)
  plot(r_slp_f)
  
  # Map of indicators for risk distribution exceeding thresholds
  r_slp_ind1 <- calc_sig_map(risk_slp_f, threshold = 0.064)
  r_slp_ind2 <- calc_sig_map(risk_slp_f, threshold = 0.032)
  r_slp_ind <- stack(r_slp_ind1, r_slp_ind2)
  plot(r_slp_ind)

  # Save rasters
  saveRDS(r_slp_f, SLP_FINE_RASTER_FILE)
  saveRDS(r_slp_ind, SLP_SIGNIFICANCE_RASTER_FILE)
  
  # Parameter estimates
  slp_params <- c()
  for (i in 1:3){
    slp_params <- rbind(
      slp_params,
      data.frame(
        Parameter = paste0('Beta (+,', i, ')'),
        Posterior.Mean = round(unname(mean(beta_pos_samples[,i])), 3),
        Posterior.80.lb = round(unname(quantile(beta_pos_samples[,i], 0.025)), 3),
        Posterior.80.ub = round(unname(quantile(beta_pos_samples[,i], 0.975)), 3)
      )
    )
  }
  for (i in 1:3){
    slp_params <- rbind(
      slp_params,
      data.frame(
        Parameter = paste0('Beta (-,', i, ')'),
        Posterior.Mean = round(unname(mean(beta_neg_samples[,i])), 3),
        Posterior.80.lb = round(unname(quantile(beta_neg_samples[,i], 0.025)), 3),
        Posterior.80.ub = round(unname(quantile(beta_neg_samples[,i], 0.975)), 3)
      )
    )
  }
  for (i in 1:3){
    slp_params <- rbind(
      slp_params,
      data.frame(
        Parameter = paste0('Beta loc (', i, ')'),
        Posterior.Mean = round(unname(mean(beta_loc_samples[,i])), 3),
        Posterior.80.lb = round(unname(quantile(beta_loc_samples[,i], 0.025)), 3),
        Posterior.80.ub = round(unname(quantile(beta_loc_samples[,i], 0.975)), 3)
      )
    )
  }
  slp_params <- rbind(
    slp_params,
    data.frame(
      Parameter = paste0('Alpha (+)'),
      Posterior.Mean = round(unname(mean(alpha_pos_samples)), 3),
      Posterior.80.lb = round(unname(quantile(alpha_pos_samples, 0.025)), 3),
      Posterior.80.ub = round(unname(quantile(alpha_pos_samples, 0.975)), 3)
    ),
    data.frame(
      Parameter = paste0('Alpha (-)'),
      Posterior.Mean = round(unname(mean(alpha_neg_samples)), 3),
      Posterior.80.lb = round(unname(quantile(alpha_neg_samples, 0.025)), 3),
      Posterior.80.ub = round(unname(quantile(alpha_neg_samples, 0.975)), 3)
    ))
  
  write_latex_table(slp_params, '/slp_estimates.txt', getwd())
  
}
