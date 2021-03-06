## Fit statistical model to smooth the raw plot biomass via cross-validation
## to determine best upper-bound on amount of spatial smoothing.

## The model in principle fits in two parts - first the proportion of points occupied by trees
## then the average biomass for occupied points (called potential biomass).
## Estimated biomass is the product of occupancy and potential.
## For total biomass, all FIA analysis done on forested plots, so occupancy is always 1.

library(dplyr)
library(assertthat)

load(file.path(interim_results_dir, paste0('full_trees_with_biomass_grid',
                                           ifelse(use_agb, '_agb', ''), '.Rda')))

library(doParallel)
if(n_cores == 0) {
    if(Sys.getenv("SLURM_JOB_ID") != "") {
        n_cores <- Sys.getenv("SLURM_CPUS_PER_TASK")
    } else n_cores <- detectCores()
}
registerDoParallel(cores = n_cores)

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

## Set up cross-validation folds.
set.seed(1)
cells <- sample(unique(cell_full$cell), replace = FALSE)
folds <- rep(1:n_folds, length.out = length(cells))

cell_full <- cell_full %>% inner_join(data.frame(cell = cells, fold = folds), by = c('cell'))

pred_occ <- matrix(1, nrow(cell_full), 1)

pred_pot_arith <- pred_pot_larith <- sig2_arith <- sig2_larith <-
    matrix(0, nrow(cell_full), length(k_pot_cv))
dimnames(pred_pot_arith)[[2]] <- dimnames(pred_pot_larith)[[2]] <-
    dimnames(sig2_arith)[[2]] <- dimnames(sig2_larith)[[2]] <- k_pot_cv

draws_logpot_arith <- draws_logpot_larith <- array(0, c(nrow(cell_full), length(k_pot_cv), n_stat_samples))

## Fit statistical model to each fold.
n_folds <- max(cell_full$fold)
output <- foreach(i = seq_len(n_folds)) %dopar% {
    train <- cell_full %>% filter(fold != i)
    test <- cell_full %>% filter(fold == i)
    po <- NULL
    ppa <- fit(train, newdata = test, k_pot = k_pot_cv, type_pot = 'arith', unc = TRUE, use_bam = TRUE, save_draws = TRUE, num_draws = n_stat_samples)
    ppl <- fit(train, newdata = test, k_pot = k_pot_cv, type_pot = 'log_arith', unc = TRUE, use_bam = TRUE, save_draws = TRUE, num_draws = n_stat_samples)
    cat("n_fold: ", i, " ", date(), "\n")
    list(po = po, ppa = ppa, ppl = ppl)
}
for(i in seq_len(n_folds)) {
    subn <- sum(cell_full$fold == i)
    pred_pot_arith[cell_full$fold == i, ] <- output[[i]]$ppa$pred_pot
    pred_pot_larith[cell_full$fold == i, ] <- output[[i]]$ppl$pred_pot
    draws_logpot_arith[cell_full$fold == i,  , ] <- output[[i]]$ppa$draws_logpot
    draws_logpot_larith[cell_full$fold == i,  , ] <- output[[i]]$ppl$draws_logpot
    ## when model not returned, we do get back the 'sig2' value as the model_pot
    sig2_arith[cell_full$fold == i, ] <- rep(output[[i]]$ppa$model_pot, each = subn)
    sig2_larith[cell_full$fold == i, ] <- rep(output[[i]]$ppl$model_pot, each = subn)
}

## Assess results

y <- cell_full$avg
## weight by points_occ as same as points_total
crit_arith <- calc_point_criterion(pred_occ, pred_pot_arith, cell_full$points_occ,
                               y, cv_max_biomass)
crit_larith <- calc_point_criterion(pred_occ, pred_pot_larith, cell_full$points_occ,
                                  y, cv_max_biomass)
colnames(crit_arith) <- colnames(crit_larith) <- k_pot_cv

cell_full$obs <- y

crit_arith <- c(list(point = crit_arith),
                     calc_cov_criterion(draws_logpot_arith, sig2 = sig2_arith,
                                        data = cell_full, type_pot = 'arith'))
crit_larith <- c(list(point = crit_larith),
                     calc_cov_criterion(draws_logpot_larith, sig2 = sig2_larith,
                                        data = cell_full, type_pot = 'log_arith'))


save(crit_arith, crit_larith, pred_pot_arith, pred_pot_larith, file = file.path(interim_results_dir, paste0('cv_total_biomass', ifelse(use_agb, '_agb', ''), '.Rda')))

