## Fit primary MCMC for a given region. Then draw composition proportions.

## Run-time: 

library(ncdf4)

region <- composition_region ## from config - don't loop over regions due to length of MCMC run

load(file.path(interim_results_dir, paste0('composition_raw_', region, '.Rda')))

latentNcdfName <- paste0('composition_fitted_latent_', region,
                         ifelse(is.null(runName), '', paste0('_', runName)), '.nc')
    
## Draw composition proportions based on latent variable model.
outputNcdfName <- paste0('composition_fitted_', region, ifelse(is.null(runName), '', paste0('_', runName)), '.nc')

make_albers_netcdf(name = 'proportion', units = 'unitless (proportion from 0 to 1)',
                   longname = 'relative composition, relative to all tree taxa,',
                   fn = outputNcdfName, dir = interim_results_dir, x = sort(unique(coord$x)),
                   y = sort(unique(coord$y)), taxa = taxa$taxonName, num_samples = floor(S/thin))

latentNcdfPtr <- nc_open(file.path(interim_results_dir, latentNcdfName))
outputNcdfPtr <- nc_open(file.path(interim_results_dir, outputNcdfName), write = TRUE)

set.seed(seed)
## this draws the proportions based on the draws of the latent variables
drawProportions(latentNcdfPtr, outputNcdfPtr, numMCsamples = numSamplesForProps,
                numInputSamples = floor(S/thin), secondThin = 1, nCells = m1*m2, taxa = taxa$taxonName)
## do draw proportions for burnin part so can assess burnin on proportions

nc_close(latentNcdfPtr)
nc_close(outputNcdfPtr)


