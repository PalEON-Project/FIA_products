## Outputs a netCDF file with dimensions of x, y, samples and
## variables for total biomass and biomass of each taxon

## need to revise for composition
## read in east and west
## blend at midpoint of OH

stop("file not set up yet")
stop("need to deal with burnin")

library(ncdf4)

# subsets to post-burnin samples

# also flips N-S so that values are increasing from S to N
# interim netCDF files go from N to S becaus that is how the tif files
# are set up

runID <- paste0(domain, "_", runID)

load(file.path(dataDir, paste0('data_', runID, '.Rda')))

outputNcdfName <- paste0('FIAcomposition_', runID, '_full.nc')

finalNcdfName <- paste0('FIA_composition_', domain, '_v', productVersion, '.nc')

outputNcdfPtr <- nc_open(file.path(outputDir, outputNcdfName))

x_grid <- sort(unique(pred_grid$x))
y_grid <- sort(unique(pred_grid$y))

x_res <- length(x_grid)
y_res <- length(y_grid)

if(domain == "western") 
    makeAlbersNetCDF(name = 'proportion', units = 'unitless (proportion from 0 to 1)',
                     longname = 'relative composition, relative to all tree taxa,',
                     fn = finalNcdfName, dir = outputDir, x = xGrid, y = yGrid,
                     taxa = taxa$taxonName, numSamples = floor((S-burnin)/(thin*secondThin)))
if(domain == "eastern")
    makeAlbersNetCDF(name = 'proportion', units = 'unitless (proportion from 0 to 1)',
                     longname = 'relative composition, relative to all tree taxa,',
                     fn = finalNcdfName, dir = outputDir, x = xGrid, y = yGrid,
                     taxa = taxa$taxonName, numSamples = floor((S-burnin)/(thin*secondThin)))

finalNcdfPtr <- nc_open(file.path(outputDir, finalNcdfName), write = TRUE)

currentSamples <- floor(S/(thin*secondThin))
finalSamples <- floor((S-burnin)/(thin*secondThin))
for(p in seq_len(nTaxa)) {
    tmp <- ncvar_get(outputNcdfPtr, varid = taxa$taxonName[p], start = c(1, 1, (currentSamples - finalSamples + 1)), count = c(-1, -1, finalSamples))
    # m2:1 flips the y-axis so that y values go from S to N
    ncvar_put(finalNcdfPtr, taxa$taxonName[p], tmp[ , m2:1, ], start = c(1, 1, 1), count = c(-1, -1, finalSamples))
}

nc_close(outputNcdfPtr)
nc_close(finalNcdfPtr)

