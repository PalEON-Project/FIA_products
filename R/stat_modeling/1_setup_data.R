library(raster)
library(dplyr)
library(readr)
library(ncdf4)


load(file.path(interim_results_dir, 'full_trees_with_biomass.Rda'))

## create base raster in albers projection
base_raster <- raster(xmn = -71000, xmx = 2297000, ncols = 296,
                    ymn = 58000,  ymx = 1498000, nrows = 180,
                    crs = '+init=epsg:3175')
base_raster <- setValues(base_raster, 1:ncell(base_raster))


if(!'cell' %in% names(fia)) 
    fia <- fia %>% add_paleon_grid(conversions_data_dir, fia_plot_to_paleon_gridcell_file, states = states, base_raster = base_raster)

## Deal with missing cell for one Rhode Island plot
if(FALSE) {  ## determines nominal cell assignment of missing RI cell
    tmp <- as.data.frame(fia_recent[obs_without_cells[1], ])
    tmp <- SpatialPointsDataFrame(coords=data.frame(x=tmp$LON, y=tmp$LAT),data=data.frame(CN=tmp$PLT_CN, STATECD=tmp$STATECD), proj4string=CRS("+init=epsg:4269"))
    tmp <- spTransform(tmp, CRS("+init=epsg:3175"))
    missing_cell <- extract(base_raster, tmp@coords)
} else missing_cell <- 32522

obs_without_cells <- which(is.na(fia$cell))
missing_loc <- unique(fia$PLT_CN[obs_without_cells])
if(length(missing_cell) != 1 || missing_loc != "122556769010661")
    stop("Unexpected missing cells") else fia$cell[obs_without_cells] <- missing_cell 

## Perhaps move FIA spcd to Paleon level 3 conversion here, as done in PLS code


xy <- xyFromCell(base_raster, getValues(base_raster))
grid <- tibble(x = xy[,1], y = xy[,2], cell = getValues(base_raster))

mask <- nc_open(file.path(conversions_data_dir, 'paleonmask.nc'))
regions <- ncvar_get(mask, 'subregion', c(1,1),c(-1,-1))

x <- matrix(ncvar_get(mask, 'x', c(1),c(-1)), nrow = nrow(regions),
            ncol = ncol(regions))
y <- matrix(ncvar_get(mask, 'y', c(1),c(-1)),nrow = nrow(regions),
            ncol = ncol(regions), byrow = TRUE) 

west <- c(regions %in% paleon_regions_west)
paleon <- c(regions %in% paleon_regions)
pred_grid <- data.frame(x = c(x), y = c(y))
pred_grid_paleon <- pred_grid[paleon, ]
pred_grid_west <- pred_grid[west, ]
