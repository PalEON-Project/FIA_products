make_albers_netcdf <- function(name = NULL, units = '', longname = '', fn = NULL, dir = '', x, y, taxa, num_samples) {
    ## Make empty netCDF files for final output
    if(is.null(fn)) stop("makeAlbersNetCDF: requires filename as argument.")
    
    x_dim <-  ncdim_def("x", "meters_east", x,
                        longname = 'x coordinate of grid cell centroid in Albers projection (Great Lakes St Lawrence Albers [Proj4 +init=epsg:3175])')
    y_dim <-  ncdim_def("y", "meters_north", y,
                        longname = 'y coordinate of grid cell centroid in Albers projection (Great Lakes St Lawrence Albers [Proj4 +init=epsg:3175])')
    dim_list <- list(x_dim, y_dim)
    if(num_samples) {
        draw_dim <- ncdim_def("sample", "number", 1:num_samples, longname = "MCMC sample")
        dim_list[[3]] <- draw_dim
    }
    vars <- list()
    length(vars) <- length(taxa)
    for(k in seq_along(taxa)) {
        vars[[k]] <- ncvar_def(name = taxa[k], dim = dim_list, units = units,
                               longname = paste0(longname, " for ", ifelse(taxa[k] == "Total", "", "taxon "),
                                                 taxa[k]), prec="double")
    }
    ncdfPtr <- nc_create(file.path(dir, fn), vars, force_v4 = TRUE)
    nc_close(ncdfPtr)
    invisible(NULL)
}
    
