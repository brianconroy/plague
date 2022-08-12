
# Preferential Sampling in Disease Surveillance

This repository contains simulations and a case study to estimate disease risk from preferentially sampled disease surveillance data monitoring plague in sciurid rodents.

## Organization

## case_study

Fits and evaluates various models to a surveillance dataset monitoring the incidence of plague in sciurid rodents throughout the state of California.

* fit_poisson.R: fits a Poisson regression model.
* fit_spatial_poisson.R: fits a spatial Poisson regression model with Gausssian process residuals.
* fit_slp.R: fits a shared latent process model attempting to correct for preferential sampling.
* cdph_results.R: summarizes results from model fitting.
* sim_from_fitted.R: simulates data from the fitted SLP model.

## simulation

* sim_datasets.R: simulates preferentially sampled disease surveillance datasets.
* sim_results.R: evaluates model performance on simulated data.

## stan

* poisson.stan: STAN file defining the Poisson model.
* spatial_poisson.stan: STAN file defining the spatial Poisson model.
* slp.stan: STAN file defining the shared latent process model.
