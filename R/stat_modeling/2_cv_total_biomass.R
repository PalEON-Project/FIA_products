## Fit statistical model to smooth the raw plot biomass.
## The model fits in two parts - first the proportion of points occupied by trees
## (this is much more important for the taxon-level fitting)
## then the average biomass for occupied points (called potential biomass).
## Estimated biomass is the product of occupancy and potential.

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

if(nrow(cell_full) != nrow(plots_per_cell))
    stop("number of cells with FIA plots not the same as the number of cells with biomass")

set.seed(1)
cells <- sample(unique(cell_full$cell), replace = FALSE)
folds <- rep(1:n_folds, length.out = length(cells))

cell_full <- cell_full %>% inner_join(data.frame(cell = cells, fold = folds), by = c('cell'))

if(n_cores == 0) {
    if(Sys.getenv("SLURM_JOB_ID") != "") {
        n_cores <- Sys.getenv("SLURM_CPUS_PER_TASK")
    } else n_cores <- detectCores()
}

library(doParallel)
registerDoParallel(cores = n_cores)

results <- fit_cv_total(cell_full, k_pot = k_pot_cv)

## assess results
    
y <- cell_full$avg

## weight by points_occ as same as points_total
critArith <- calc_cv_criterion(results$pred_occ, results$pred_pot_arith, cell_full$points_occ, y, cv_max_biomass)
critLogArith <- calc_cv_criterion(results$pred_occ, results$pred_pot_larith, cell_full$points_occ, y, cv_max_biomass)
colnames(critArith) <- colnames(critLogArith) <- k_pot_cv

save(critArith, critLogArith, results, file = file.path(interim_results_dir, 'cv_total_biomass.Rda'))

