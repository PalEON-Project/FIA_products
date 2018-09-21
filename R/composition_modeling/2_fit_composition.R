library(ncdf4)

region <- composition_region ## from config - don't loop over regions due to length of MCMC run

if(!secondRun) {
    latentNcdfName <- paste0('composition_fitted_latent_', region, '.nc')
} else latentNcdfName <- paste0('composition_fitted_latent_run2_', region, '.nc')
    
if(do_mcmc) {
    load(file.path(interim_results_dir, paste0('composition_raw_', region, '.Rda')))
    
    if(!resumeRun && !secondRun) 
        set.seed(seed)

    if(!resumeRun)
        make_albers_netcdf(name = 'latent', units = 'unitless', longname = 'latent multivariate logit values', fn = latentNcdfName, dir = interim_results_dir, x = sort(unique(coord$x)), y = sort(unique(coord$y)), taxa = taxa$taxonName, num_samples = floor(S/thin))
    
    out <- runMCMC(y = data$taxon, cell = data$cell, plot = data$plot, plotInfo = plotInfo, gridInfo = coord, C = nbhd, Cindices = nbhdIndices, S = S, thin = thin, resumeRun = resumeRun, secondRun = secondRun, hyperpar = c(-0.5, 0), nbhdStructure = nbhdStructure, outputNcdfName = latentNcdfName, taxa = taxa, adaptStartDecay = ifelse(secondRun, 0, burnin), runID = paste0(region, '-', runID, sep = ''), dataDir = interim_results_dir, outputDir = interim_results_dir)
}

if(do_draws) {
    outputNcdfName <- paste0('composition_fitted_', region, '.nc')

    make_albers_netcdf(name = 'proportion', units = 'unitless (proportion from 0 to 1)', longname = 'relative composition, relative to all tree taxa,', fn = outputNcdfName, dir = interim_results_dir, x = sort(unique(coord$x)), y = sort(unique(coord$y)), taxa = taxa$taxonName, num_samples = floor(S/thin))

    latentNcdfPtr <- nc_open(file.path(interim_results_dir, latentNcdfName))
outputNcdfPtr <- nc_open(file.path(interim_results_dir, outputNcdfName), write = TRUE)

    ## this draws the proportions based on the draws of the latent variables
    set.seed(seed)

    drawProportions(latentNcdfPtr, outputNcdfPtr, numMCsamples = numSamplesForProps, numInputSamples = floor(S/thin), secondThin = 1, nCells = m1*m2, taxa = taxa$taxonName)
    ## do draw proportions for burnin part so can assess burnin on proportions

    nc_close(latentNcdfPtr)
    nc_close(outputNcdfPtr)
}
