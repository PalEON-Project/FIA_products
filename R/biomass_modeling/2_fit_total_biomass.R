## Fit statistical model to smooth the raw plot biomass.
## The model in principle fits in two parts - first the proportion of points occupied by trees
## then the average biomass for occupied points (called potential biomass).
## Estimated biomass is the product of occupancy and potential.
## For total biomass, all FIA analysis done on forested plots, so occupancy is always 1.

load(file.path(interim_results_dir, 'full_trees_with_biomass_grid.Rda'))

if(!exists('k_pot_total'))
    stop("Must specify 'k_pot_total'")

## for naming consistency with PLS, will refer to FIA plots as 'points'


area_conversion <- acres_per_ha * squ_feet_per_acre /
    (pi * fia_subplot_radius_in_feet^2 * fia_subplots_per_plot) / kg_per_Mg 

## Total biomass per plot.
biomass_plot <- fia %>% group_by(PLT_CN, x, y, cell) %>%
    summarize(total_biomass = sum(biomass) * area_conversion)
## after conversion, biomass is Mg/Ha

## Average biomass per cell.
## Note that for total biomass, all plots are 'occupied' as min biomass > 0
## essentially by definition since we exclude non-forested plots
cell_full <- biomass_plot %>% group_by(cell, x, y) %>%
    summarize(points_occ = n(),
              avg = mean(total_biomass),
              geom_avg = mean(log(total_biomass))
              )

## Total number of fia points per grid cell
plots_per_cell <- biomass_plot %>% group_by(cell) %>% summarize(points_total = n()) 

assert_that(nrow(cell_full) == nrow(plots_per_cell),
            msg = "Number of cells with FIA plots not the same as the number of cells with biomass")

## fit stats model, using potential biomass part of model only.
biomass_total <- fit(cell_full, newdata = pred_grid_paleon, k_occ = NULL,
                     k_pot = k_pot_total, return_model = TRUE, unc = TRUE,
                     type_pot = 'arith', num_draws = n_stat_samples, save_draws = TRUE,
                     use_bam = TRUE)

save(biomass_total, file = file.path(interim_results_dir, 'fitted_total_biomass.Rda'))

