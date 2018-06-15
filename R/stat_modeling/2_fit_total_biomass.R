## Fit statistical model to smooth the raw plot biomass.
## The model fits in two parts - first the proportion of points occupied by trees
## (this is much more important for the taxon-level fitting)
## then the average biomass for occupied points (called potential biomass).
## Estimated biomass is the product of occupancy and potential.

## Run-time (without cross-validation) approximately 5 minutes with k value of 1000.

if(!exists('k_pot_total'))
    stop("Must specify 'k_pot_total'")



## for naming consistency with PLS, will refer to FIA plots as 'points'

area_conversion <- acres_per_ha * squ_feet_per_acre /
    (pi * fia_subplot_radius_in_feet^2 * fia_subplots_per_plot) / kg_per_Mg 

biomass_plot <- fia %>% group_by(PLT_CN, x, y, cell) %>%
    summarize(total_biomass = sum(biomass) * area_conversion)
## after conversion, biomass is Mg/Ha

## note that for total biomass, all plots are 'occupied' as min biomass > 0
## essentially by definition since we exclude non-forested plots
cell_full <- biomass_plot %>% group_by(cell, x, y) %>%
    summarize(points_occ = n(),
              avg = mean(total_biomass),
              geom_avg = mean(log(total_biomass))
              )

## total number of fia points per grid cell
plots_per_cell <- biomass_plot %>% group_by(cell) %>% summarize(points_total = n()) 

## TODO: check for natural grasslands/plots without tree
## and whether such plots would make it through the filters in the code
if(nrow(cell_full) != nrow(plots_per_cell))
    stop("number of cells with FIA plots not the same as the number of cells with biomass")

if(do_cv) {
    stop("code for CV not fully finished")
    k_pot = c(100,250,500,1000,1500,2000,2500,3000,3500)
    
    set.seed(1)
    cells <- sample(unique(cell_full$cell), replace = FALSE)
    folds <- rep(1:n_folds, length.out = length(cells))
    
    cell_full <- cell_full %>% inner_join(data.frame(cell = cells, fold = folds), by = c('cell'))
    stop("this CV code won't work yet")
    results <- fit_cv(cell_full, k_occ, k_pot, n_cores)

    ## assess results
    
    y <- cell_full$avg*cell_full$points_occ/cell_full$points_total ## actual average biomass over all cells
    y[is.na(y)] <- 0
    y[y > mx] <- mx


    critArith <- calc_cv_criterion(results$pred_occ, results$pred_pot_arith, cell_full$points_total,
                                   cell_full$avg*cell_full$count/cell_full$total, 200)
    critLogArith <- calc_cv_criterion(results$pred_occ, results$pred_pot_larith, cell_full$points_total,
                                   cell_full$avg*cell_full$count/cell_full$total, 200)

}

## fit stats model
biomass_total <- fit(cell_full, newdata = pred_grid_west, k_occ = NULL, k_pot = k_pot_total, return_model = TRUE, unc = TRUE, type_pot = 'log_arith', num_draws = n_stat_samples)

save(biomass_total, file = file.path(interim_results_dir, 'fitted_total_biomass.Rda'))

