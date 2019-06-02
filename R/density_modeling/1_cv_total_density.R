## Fit statistical model to smooth the raw plot density via cross-validation
## to determine best upper-bound on amount of spatial smoothing.

## The model in principle fits in two parts - first the proportion of points occupied by trees
## then the average density for occupied points (called potential density).
## Estimated density is the product of occupancy and potential.
## For total density, all FIA analysis done on forested plots, so occupancy is always 1.

library(dplyr)

load(file.path(interim_results_dir, 'full_trees_with_biomass_grid.Rda'))

library(doParallel)
if(n_cores == 0) {
    if(Sys.getenv("SLURM_JOB_ID") != "") {
        n_cores <- Sys.getenv("SLURM_CPUS_PER_TASK")
    } else n_cores <- detectCores()
}
registerDoParallel(cores = n_cores)

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

## Set up cross-validation folds.
set.seed(1)
cells <- sample(unique(cell_full$cell), replace = FALSE)
folds <- rep(1:n_folds, length.out = length(cells))

cell_full <- cell_full %>% inner_join(data.frame(cell = cells, fold = folds), by = c('cell'))

pred_occ <- matrix(1, nrow(cell_full), 1)

pred_pot_arith <- pred_pot_larith <- matrix(0, nrow(cell_full), length(k_pot_cv))
dimnames(pred_pot_arith)[[2]] <- dimnames(pred_pot_larith)[[2]] <- k_pot_cv
    
## Fit statistical model to each fold.
n_folds <- max(cell_full$fold)
output <- foreach(i = seq_len(n_folds)) %dopar% {
    train <- cell_full %>% filter(fold != i)
    test <- cell_full %>% filter(fold == i)
    po <- NULL
    ppa <- fit(train, newdata = test, k_pot = k_pot_cv, type_pot = 'arith', unc = FALSE, use_bam = TRUE)
    ppl <- fit(train, newdata = test, k_pot = k_pot_cv, type_pot = 'log_arith', unc = FALSE, use_bam = TRUE)
    cat("n_fold: ", i, " ", date(), "\n")
    list(po, ppa, ppl)
}
for(i in seq_len(n_folds)) {
    pred_pot_arith[cell_full$fold == i, ] <- output[[i]][[2]]$pred_pot
    pred_pot_larith[cell_full$fold == i, ] <- output[[i]][[3]]$pred_pot
}

## Assess results

y <- cell_full$avg
## weight by points_occ as same as points_total
critArith <- calc_point_criterion(pred_occ, pred_pot_arith, cell_full$points_occ,
                               y, cv_max_density)
critLogArith <- calc_point_criterion(pred_occ, pred_pot_larith, cell_full$points_occ,
                                  y, cv_max_density)
colnames(critArith) <- colnames(critLogArith) <- k_pot_cv

save(critArith, critLogArith, pred_pot_arith, pred_pot_larith, file = file.path(interim_results_dir, 'cv_total_density.Rda'))

