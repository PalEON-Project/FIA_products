## Fit primary MCMC for a given region. Then draw composition proportions.

library(ncdf4)

region <- composition_region ## from config - don't loop over regions due to length of MCMC run

load(file.path(interim_results_dir, paste0('composition_raw_', region, '.Rda')))


if(secondRun && resumeRun)
    stop("Not clear if code is properly set up to resume from middle of second run.")

if(!secondRun) {
    latentNcdfName <- paste0('composition_fitted_latent_', region, '.nc')
} else latentNcdfName <- paste0('composition_fitted_latent_run2_', region, '.nc')
    
load(file.path(interim_results_dir, paste0('composition_raw_', region, '.Rda')))

if(!resumeRun && !secondRun) 
    set.seed(seed)

if(!resumeRun)
    make_albers_netcdf(name = 'latent', units = 'unitless', longname = 'latent multivariate logit values',
                       fn = latentNcdfName, dir = interim_results_dir,
                       x = sort(unique(coord$x)), y = sort(unique(coord$y)),
                       taxa = taxa$taxonName, num_samples = floor(S/thin))

out <- runMCMC(y = data$taxon, cell = data$cell, plot = data$plot, plotInfo = plotInfo,
               gridInfo = coord, C = nbhd, Cindices = nbhdIndices, S = S, thin = thin,
               resumeRun = resumeRun, secondRun = secondRun, hyperpar = c(-0.5, 0),
               nbhdStructure = nbhdStructure, outputNcdfName = latentNcdfName, taxa = taxa,
               adaptStartDecay = ifelse(secondRun, 0, burnin),
               runID = paste0(region, '-', runID, sep = ''),
               dataDir = interim_results_dir, outputDir = interim_results_dir)


