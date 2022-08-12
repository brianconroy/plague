
# Data Dictionary

Due to data sharing restrictions from the California Department of Public Health
stemming from privacy concerns around the implications of plague occurrence for 
real estate prices, etc, the counts of positive and negative rodents in this dataset
have been simulated from our fitted model. The R data file accompanying this paper contains a list with the following elements.

* **Y_loc**: counts of sampling events conducted in each grid cell of the study region.
* **Y_pos**: counts of disease positive rodents collected in each grid cell that has been sampled by the surveillance system.
* **Y_neg**: counts of negative positive rodents collected in each grid cell that has been sampled by the surveillance system.
* **X**: matrix of covariates and intercept along with the first and second principal components of the PRISM climatic dataset.
* **prism_rasters**: rasters of the first and second principal components of the PRISM climatic dataset.

The ordering of values in Y_loc, Y_pos, and Y_neg corresponds to the indices of grid cells in the PRISM rasters. For instance, if you wish to visualize counts of sampling events, or of collelcted rodents, on the raster defining the study region, apply the following:

```
library(raster)
library(sp)

plague_data <- readRDS('plague_data_simulated.rds')
prism_rasters <- plague_data$prism_rasters
Y_loc <- plague_data$Y_loc
Y_pos <- plague_data$Y_pos

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

```