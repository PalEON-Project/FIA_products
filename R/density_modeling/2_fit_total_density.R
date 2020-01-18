## Fit statistical model to smooth the raw plot density.
## The model in principle fits in two parts - first the proportion of points occupied by trees
## then the average density for occupied points (called potential density).
## Estimated density is the product of occupancy and potential.
## For total density, all FIA analysis done on forested plots, so occupancy is always 1.

library(dplyr)
library(assertthat)

load(file.path(interim_results_dir, 'full_trees_with_biomass_grid.Rda'))

if(!exists('k_pot_total_density'))
    stop("Must specify 'k_pot_total_density'")

## for naming consistency with PLS, will refer to FIA plots as 'points'


area_conversion <- acres_per_ha * squ_feet_per_acre /
    (pi * fia_subplot_radius_in_feet^2 * fia_subplots_per_plot)

## Total density per plot.
density_plot <- fia %>% group_by(PLT_CN, x, y, cell) %>%
    summarize(total_density = n() * area_conversion)
## after conversion, density is stems/Ha

## Average density per cell.
## Note that for total density, all plots are 'occupied' as min density > 0
## essentially by definition since we exclude non-forested plots
cell_full <- density_plot %>% group_by(cell, x, y) %>%
    summarize(points_occ = n(),
              avg = mean(total_density),
              geom_avg = mean(log(total_density))
              )

## Total number of fia points per grid cell
plots_per_cell <- density_plot %>% group_by(cell) %>% summarize(points_total = n()) 

assert_that(nrow(cell_full) == nrow(plots_per_cell),
            msg = "Number of cells with FIA plots not the same as the number of cells with density")

## fit stats model, using potential density part of model only.
density_total <- fit(cell_full, newdata = pred_grid_paleon, k_occ = NULL,
                     k_pot = k_pot_total_density, return_model = TRUE, unc = TRUE,
                     type_pot = fit_scale_density, num_draws = n_stat_samples, save_draws = TRUE,
                     use_bam = TRUE)

save(density_total, file = file.path(interim_results_dir, 'fitted_total_density.Rda'))

