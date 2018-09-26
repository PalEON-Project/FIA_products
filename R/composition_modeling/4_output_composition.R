## Outputs a netCDF file with dimensions of x, y, samples and
## variables for total biomass and biomass of each taxon

## need to revise for composition
## read in east and west
## blend at midpoint of OH

library(ncdf4)

for(region in c('eastern', 'western')) {
    load(file.path(interim_results_dir, paste0('composition_raw_', region, '.Rda')))
    assign(paste0('taxa_', region), taxa)
    output_ncdf_name <- paste0('composition_fitted_', region,
                               ifelse(first_run_name != '', paste0('_', first_run_name), ''), '.nc')
    if(second_run_name != '')
        output_ncdf_name <- c(output_ncdf_name,
                            paste0('composition_fitted_', region, '_', second_run_name, '.nc'))
    assign(paste0('output_ncdf_name_', region), output_ncdf_name)
}

assert_that(identical(taxa_eastern, taxa_western), msg = 'mismatched taxa in regions')
taxa <- taxa_eastern


mask <- nc_open(file.path(conversions_data_dir, 'paleonmask.nc'))
regions <- ncvar_get(mask, 'subregion', c(1,1), c(-1,-1))
domain <- ncvar_get(mask, 'domain', c(1,1), c(-1,-1))

x <- matrix(ncvar_get(mask, 'x', c(1),c(-1)), nrow = nrow(regions),
            ncol = ncol(regions))
y <- matrix(ncvar_get(mask, 'y', c(1),c(-1)),nrow = nrow(regions),
            ncol = ncol(regions), byrow = TRUE) 

## full rectangular grid
pred_grid <- data.frame(x = c(x), y = c(y))
## subset to only land areas within PalEON states:
pred_grid_paleon <- pred_grid[regions %in% paleon_regions, ]

x_grid <- sort(unique(pred_grid$x))
y_grid <- rev(sort(unique(pred_grid$y)))

x_res <- length(x_grid)
y_res <- length(y_grid)

final_ncdf_name <- paste0('FIA_composition_v', product_version, '.nc')

make_albers_netcdf(name = 'proportion', units = 'unitless (proportion from 0 to 1)',
                 longname = 'relative composition, relative to all tree taxa,',
                 fn = final_ncdf_name, dir = output_dir, x = xGrid, y = yGrid,
                 taxa = taxa$taxonName, num_samples = floor((S-burnin)/(thin*secondThin)))

final_ncdf_ptr <- nc_open(file.path(output_dir, final_ncdf_name), write = TRUE)

## The blending happens at a point in Ohio such that the vertical line does not touch
## Michigan and such that the western edge of the blending area is east of where
## Ohio and Michigan meet and the eastern edge of the blending area is west of the Ohio-
## Pennsylvania border. Also, blend_below needs to be set such that areas in Michigan
## that are east of the western edge of the blending area are excluded from blending.

## To illustrate:
if(FALSE) {
    plot(pred_grid_paleon$x, pred_grid_paleon$y, pch = 16, cex = 0.6)
    abline(h = y_grid[blend_below])
    abline(v = x_grid[blend_point - blend_buffer])
    abline(v = x_grid[blend_point + blend_buffer])
}

if(blend_point - blend_buffer < 136 || blend_point + blend_buffer > 159)
    stop("code for blending east and west not set up for current blending area")

## Create weight matrices to do blending. These are created with rows as y and columns as x
## to make it easier to think about in same orientation as a map.
## Counting of y cells goes north to south for the moment to match values in interim netCDF files
x_vals <- matrix(rep(1:x_res, each = y_res), nrow = y_res, ncol = x_res)
y_vals <- t(matrix(rep(1:y_res, each = x_res), nrow = x_res, ncol = y_res))

wgts <- seq(1, 0, length = blend_buffer*2 + 1)  # linear blending of east and west
blend_ids <- (blend_point - blend_buffer):(blend_point + blend_buffer)
num_below <- sum(y_vals[ , 1] > blend_below)

wgts_west <- wgts_east <- matrix(0, y_res, x_res)
wgts_west[y_vals <= blend_below & x_vals < blend_point] <- 1
wgts_west[y_vals > blend_below & x_vals < blend_point - blend_buffer] <- 1
wgts_west[y_vals > blend_below & x_vals %in% blend_ids] <- rep(wgts, each = num_below)

wgts_east[x_vals >= blend_point] <- 1
wgts_east[y_vals > blend_below & x_vals %in% blend_ids] <- rep(1-wgts, each = num_below)

## Now transpose so x is rows to match what is returned by ncvar_get().
wgts_east <- t(wgts_east) 
wgts_west <- t(wgts_west) 

if(second_run_name != "") 
    S <- S/2

output_ncdf_ptr_eastern <- nc_open(file.path(interim_results_dir, output_ncdf_name_eastern[1]))
output_ncdf_ptr_western <- nc_open(file.path(interim_results_dir, output_ncdf_name_western[1]))

current_samples <- floor(S/(thin*secondThin))
final_samples <- floor((S-burnin)/(thin*secondThin))
output <- array(0, c(x_res, y_res, final_samples))
tmp_east <- array(0, c(x_res, y_res, final_samples))
tmp_west <- array(0, c(x_res, y_res, final_samples))

for(p in seq_len(nTaxa)) {
    vals_east <- ncvar_get(output_ncdf_ptr_eastern, varid = taxa$taxonName[p],
                           start = c(1, 1, (current_samples - final_samples + 1)),
                           count = c(-1, -1, final_samples))
    vals_west <- ncvar_get(output_ncdf_ptr_western, varid = taxa$taxonName[p],
                           start = c(1, 1, (current_samples - final_samples + 1)),
                           count = c(-1, -1, final_samples))

    tmp_east[easternDomainX, easternDomainY, ] <- vals_east
    tmp_west[westernDomainX, westernDomainY, ] <- vals_west
    for(k in seq_len(final_samples)) {
        ## Blend east and west and flip the y direction so final netCDF file has y
        ## go from south to north, as in other PalEON products.
        ## Also mask to PalEON states.
        output[ , , k] <- tmp_east[ , , k] * wgts_east + tmp_west[ , , k] * wgts_west
        output[ , , k] <- domain * output[ , y_res:1, k]
    }
    ncvar_put(final_ncdf_ptr, taxa$taxonName[p], output,
              start = c(1, 1, 1), count = c(-1, -1, final_samples))
}

nc_close(output_ncdf_ptr_eastern)
nc_close(output_ncdf_ptr_western)

if(second_run_name != "") {
    start_samples <- final_samples + 1
    output_ncdf_ptr_eastern <- nc_open(file.path(interim_results_dir, output_ncdf_name_eastern[2]))
    output_ncdf_ptr_western <- nc_open(file.path(interim_results_dir, output_ncdf_name_western[2]))
    
    current_samples <- floor(S/(thin*secondThin))
    final_samples <- floor(S/(thin*secondThin))
    output <- array(0, c(x_res, y_res, final_samples))
    tmp_east <- array(0, c(x_res, y_res, final_samples))
    tmp_west <- array(0, c(x_res, y_res, final_samples))
    
    for(p in seq_len(nTaxa)) {
        vals_east <- ncvar_get(output_ncdf_ptr_eastern, varid = taxa$taxonName[p],
                               start = c(1, 1, (current_samples - final_samples + 1)),
                               count = c(-1, -1, final_samples))
        vals_west <- ncvar_get(output_ncdf_ptr_western, varid = taxa$taxonName[p],
                               start = c(1, 1, (current_samples - final_samples + 1)),
                               count = c(-1, -1, final_samples))
        
        tmp_east[easternDomainX, easternDomainY, ] <- vals_east
        tmp_west[westernDomainX, westernDomainY, ] <- vals_west
        for(k in seq_len(final_samples)) {
            ## Blend east and west and flip the y direction so final netCDF file has y
            ## go from south to north, as in other PalEON products.
            ## Also mask to PalEON states.
            output[ , , k] <- tmp_east[ , , k] * wgts_east + tmp_west[ , , k] * wgts_west
            output[ , , k] <- output[ , y_res:1, k]
            output[ , , k][domain == 0] <- NA
        }
        ncvar_put(final_ncdf_ptr, taxa$taxonName[p], output, start = c(1, 1, start_samples), count = c(-1, -1, final_samples))
    }
    
    nc_close(output_ncdf_ptr_eastern)
    nc_close(output_ncdf_ptr_western)
}

nc_close(final_ncdf_ptr)

