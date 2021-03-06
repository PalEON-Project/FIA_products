## Fit statistical model to smooth the raw plot density via cross-validation
## to determine best upper-bound on amount of spatial smoothing.

## The model fits in two parts - first the proportion of points occupied by trees
## (this is much more important for the taxon-level fitting)
## then the average density for occupied points (called potential density).
## Estimated density is the product of occupancy and potential.

## This is very computationally-intensive and best done on a cluster.

library(dplyr)

load(file.path(interim_results_dir, 'full_trees_with_biomass_grid.Rda'))

if(use_mpi) {
    library(doMPI)
    cl <- startMPIcluster()
    registerDoMPI(cl)
} else {
    library(doParallel)
    if(n_cores == 0) {
        if(Sys.getenv("SLURM_JOB_ID") != "") {
            n_cores <- Sys.getenv("SLURM_CPUS_PER_TASK")
        } else n_cores <- detectCores()
    }
    registerDoParallel(cores = n_cores)
}

taxa_to_fit <- unique(fia$level3s)
taxa_to_fit <- taxa_to_fit[!is.na(taxa_to_fit)]  ## 22 Douglas fir trees in dataset. Not fit.
print(taxa_to_fit)

## for naming consistency with PLS, will refer to FIA plots as 'points'

area_conversion <- acres_per_ha * squ_feet_per_acre /
    (pi * fia_subplot_radius_in_feet^2 * fia_subplots_per_plot)

## total number of fia points per grid cell
plots_per_cell <- fia %>% group_by(PLT_CN, x, y, cell) %>%
    summarize(n_trees = n()) %>% 
    group_by(cell) %>% summarize(points_total = n()) 

## Total density per taxon per plot only for taxa represented in a plot.
## Currently grouping taxa by level3s but this could change.
density_taxon_plot <- fia %>% group_by(PLT_CN, level3s, cell, x, y) %>%
    summarize(total_density = n() * area_conversion)
## after conversion, density is stems/Ha

## Average density per taxon per cell only for taxa represented in a cell,
## and only for plots in which taxon is present, so that potential density
## modeling ignores occupancy.
density_taxon_cell <- density_taxon_plot %>% group_by(cell, x, y, level3s) %>%
    summarize(avg = mean(total_density),
              geom_avg = mean(log(total_density)), points_occ = n()) 

## Flesh out data so have 0s for taxa missing in each cell
## (needed for occupancy part of model)
tmp <- fia %>% dplyr::select(cell, x, y) %>% unique()
expanded <- expand.grid(level3s = sort(unique(fia$level3s)),
                        cell = tmp$cell, 
                        stringsAsFactors = FALSE) %>%
            inner_join(tmp, by = 'cell')
cell_full <- expanded %>% 
    left_join(density_taxon_cell, by = c('cell', 'level3s', 'x', 'y')) %>%
    left_join(plots_per_cell, by = 'cell') %>%
    mutate(points_occ = ifelse(is.na(points_occ), 0, points_occ),
           avg = ifelse(is.na(avg), 0, avg),
           geom_avg = ifelse(is.na(geom_avg), 0, geom_avg))

assert_that(nrow(cell_full) == nrow(plots_per_cell) * length(taxa_to_fit),
            msg = "Dataset does not have data for all taxa for all occupied cells.")

## Set up cross-validation folds.
set.seed(1)
cells <- sample(unique(cell_full$cell), replace = FALSE)
folds <- rep(1:n_folds, length.out = length(cells))

cell_full <- cell_full %>% inner_join(data.frame(cell = cells, fold = folds), by = c('cell'))

## Fit statistical model to each taxon and fold.
## Nested foreach will run separate tasks for each combination of taxon and fold.
output <- foreach(taxonIdx = seq_along(taxa_to_fit)) %:%
    foreach(i = seq_len(n_folds)) %dopar% {

        taxon <- taxa_to_fit[taxonIdx]
        sub <- cell_full %>% filter(level3s == taxon)

        train <- sub %>% filter(fold != i)
        test <- sub %>% filter(fold == i)

        po <- fit(train, newdata = test, k_occ = k_occ_cv, unc = FALSE, use_bam = TRUE)
        ppa <- fit(train, newdata = test, k_pot = k_pot_cv, type_pot = 'arith', unc = FALSE, use_bam = TRUE)
        ppl <- fit(train, newdata = test, k_pot = k_pot_cv, type_pot = 'log_arith', unc = FALSE, use_bam = TRUE)
        cat("taxon: ", taxonIdx, " ; fold: ", i, "\n", sep = "")
        list(po$pred_occ, ppa$pred_pot, ppl$pred_pot)
    }

## Gather results together

pred_occ <- array(0, c(length(taxa_to_fit), nrow(plots_per_cell), length(k_occ_cv)))
pred_pot_arith <- pred_pot_larith <- array(0, c(length(taxa_to_fit), nrow(plots_per_cell), length(k_pot_cv)))

dimnames(pred_occ)[[1]] <- dimnames(pred_pot_arith)[[1]] <- dimnames(pred_pot_larith)[[1]] <- taxa_to_fit
dimnames(pred_occ)[[3]] <- k_occ_cv
dimnames(pred_pot_arith)[[3]] <- dimnames(pred_pot_larith)[[3]] <- k_pot_cv

for(taxonIdx in seq_along(taxa_to_fit))
    for(i in seq_len(n_folds)) {
        taxon <- taxa_to_fit[taxonIdx]
        sub <- cell_full %>% filter(level3s == taxon)
        pred_occ[taxonIdx, sub$fold == i, ] <- output[[taxonIdx]][[i]][[1]]
        pred_pot_arith[taxonIdx, sub$fold == i, ] <- output[[taxonIdx]][[i]][[2]]
        pred_pot_larith[taxonIdx, sub$fold == i, ] <- output[[taxonIdx]][[i]][[3]]
    }


## Assess results

critArith <- critLogArith <- array(0, c(length(taxa_to_fit), length(k_occ_cv), length(k_pot_cv)))
dimnames(critArith)[[1]] <- dimnames(critLogArith)[[1]] <- taxa_to_fit
dimnames(critArith)[[2]] <- dimnames(critLogArith)[[2]] <- k_occ_cv
dimnames(critArith)[[3]] <- dimnames(critLogArith)[[3]] <- k_pot_cv

for(taxonIdx in seq_along(taxa_to_fit)) {
    ## extract raw data (again) for the taxon
    taxon <- taxa_to_fit[taxonIdx]
    sub <- cell_full %>% filter(level3s == taxon)
        
    y <- sub$avg*sub$points_occ/sub$points_total  ## actual average density over all plots (occupied or not) in a cell

    critArith[taxonIdx, , ] <- calc_point_criterion(pred_occ[taxonIdx, , ], pred_pot_arith[taxonIdx, , ],
                                                 sub$points_total, y, cv_max_density)
    critLogArith[taxonIdx, , ] <- calc_point_criterion(pred_occ[taxonIdx, , ], pred_pot_larith[taxonIdx, , ],
                                                    sub$points_total, y, cv_max_density)
}

save(critArith, critLogArith, pred_occ, pred_pot_arith, pred_pot_larith, file = file.path(interim_results_dir, 'cv_taxon_density.Rda'))

if(use_mpi) closeCluster(cl)
