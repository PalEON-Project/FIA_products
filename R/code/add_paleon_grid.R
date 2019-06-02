library(readr)
library(dplyr)
library(sp)
library(readxl)

## Add PalEON cell info to dataset
add_paleon_grid <- function(data, conversions_data_dir, fia_plot_to_paleon_gridcell_file,
                            states, base_raster) {
    file <- file.path(conversions_data_dir, fia_plot_to_paleon_gridcell_file)
    if(!file.exists(file)) {
        grid <- add_grid_to_fia(file, base_raster)
    } else
        grid <- read_csv(file.path(conversions_data_dir, fia_plot_to_paleon_gridcell_file), 
                         col_types = 'dcddd')  %>%  ## dplyr::select(-LON_ALBERS, -LAT_ALBERS) %>%
                rename(x = cell_x, y = cell_y)
    if(!'PLT_CN' %in% names(data))
        stop("Column name of plot census ID (PLT_CN) has changed")
    if(!'CN' %in% names(grid) || !'cell' %in% names(grid))
        stop("Column names in fia to grid conversion file have changed")
    grid <- grid %>% dplyr::select(-STATECD)  # -x, -y
    return(data %>% left_join(grid, by = c("PLT_CN" = "CN")))
}

## create FIA_raster_cell_albers.csv that maps FIA plots to PalEON cells
## this uses Excel files from USFS that provide unfuzzed and unswapped Lat/Lon values
add_grid_to_fia <- function(file, base_raster) {
    stop("FIA exact location data files no longer available.")  # As of July 2019
    x1 <- read_excel(file.path(conversions_data_dir, nrs_coordinates_file)) 
    x2 <- read_excel(file.path(conversions_data_dir, nrs_coordinates_file_run2)) %>%
        select(STATECD, CN, LAT_ACTUAL_NAD83, LON_ACTUAL_NAD83)
    x3 <- read_excel(file.path(conversions_data_dir, nrs_coordinates_file_run3)) %>%
        select(STATECD, CN, LAT_ACTUAL_NAD83, LON_ACTUAL_NAD83)

    x <- rbind(x1, x2, x3)

    x <- SpatialPointsDataFrame(coords = data.frame(x = x$LON_ACTUAL_NAD83, y = x$LAT_ACTUAL_NAD83),
                                data = data.frame(CN = x$CN, STATECD = x$STATECD),
                                proj4string = CRS("+init=epsg:4269"))
    ## Transform to PalEON Albers projection
    x <-spTransform(x, CRS("+init=epsg:3175"))

    grid <- x@data
    grid$cell <- extract(base_raster, x@coords)

    write_csv(grid, path = file)
  return(grid)
}
