## Fit statistical model to smooth the raw plot biomass.
## The model fits in two parts - first the proportion of points occupied by trees
## (this is much more important for the taxon-level fitting)
## then the average biomass for occupied points (called potential biomass).
## Estimated biomass is the product of occupancy and potential.

l3a_to_l3s <- read_csv(file.path(conversions_data_dir, level3a_to_level3s_file)) 
taxa_to_fit <- unique(l3a_to_l3s$level3s)
taxa_to_fit <- taxa_to_fit[!taxa_to_fit %in% excluded_level3s_OH]
print(taxa_to_fit)

if(use_mpi) {
    library(doMPI)
    cl <- startMPIcluster()
    registerDoMPI(cl)
} else {
    if(n_cores == 0) {
        if(Sys.getenv("SLURM_JOB_ID") != "") {
            n_cores <- Sys.getenv("SLURM_CPUS_PER_TASK")
        } else n_cores <- detectCores()
    }
    library(doParallel)
    registerDoParallel(cores = n_cores)
}

## for naming consistency with PLS, will refer to FIA plots as 'points'

area_conversion <- acres_per_ha * squ_feet_per_acre /
    (pi * fia_subplot_radius_in_feet^2 * fia_subplots_per_plot) / kg_per_Mg 

## total number of fia points per grid cell
plots_per_cell <- fia %>% group_by(PLT_CN, x, y, cell) %>%
    summarize(n_trees = n()) %>% 
    group_by(cell) %>% summarize(points_total = n()) 

## currently grouping taxa by level3s but this is likely to change
biomass_taxon_plot <- fia %>% group_by(PLT_CN, level3s, cell, x, y) %>%
    summarize(total_biomass = sum(biomass) * area_conversion)
## after conversion, biomass is Mg/Ha

## here the count will be less than the number of FIA plots per cell
## when a taxon is missing from a plot
biomass_taxon_cell <- biomass_taxon_plot %>% group_by(cell, x, y, level3s) %>%
    summarize(avg = mean(total_biomass),
              geom_avg = mean(log(total_biomass)), points_occ = n()) 

## flesh out data so have 0s for taxa missing in each cell (needed for occupancy part of model)
tmp <- fia %>% dplyr::select(cell, x, y) %>% unique()
expanded <- expand.grid(level3s = sort(unique(fia$level3s)),
                        cell = tmp$cell, 
                        stringsAsFactors = FALSE) %>%
            inner_join(tmp, by = 'cell')
cell_full <- expanded %>% 
    left_join(biomass_taxon_cell, by = c('cell', 'level3s', 'x', 'y')) %>%
    left_join(plots_per_cell, by = 'cell') %>%
    mutate(points_occ = ifelse(is.na(points_occ), 0, points_occ),
           avg = ifelse(is.na(avg), 0, avg),
           geom_avg = ifelse(is.na(geom_avg), 0, geom_avg))

set.seed(1)
cells <- sample(unique(cell_full$cell), replace = FALSE)
folds <- rep(1:n_folds, length.out = length(cells))

cell_full <- cell_full %>% inner_join(data.frame(cell = cells, fold = folds), by = c('cell'))

output <- foreach(taxonIdx = seq_along(taxa_to_fit)) %:%
    foreach(i = seq_len(n_folds)) %dopar% {

        taxon <- taxa_to_fit[taxonIdx]
        sub <- cell_full %>% filter(level3s == taxon)

        train <- sub %>% filter(fold != i)
        test <- sub %>% filter(fold == i)

        po <- fit(train, newdata = test, k_occ = k_occ_cv, unc = FALSE, use_bam = TRUE)
        ppa <- fit(train, newdata = test, k_pot = k_pot_cv, type_pot = 'arith', unc = FALSE, use_bam = TRUE)
        ppl <- fit(train, newdata = test, k_pot = k_pot_cv, type_pot = 'log_arith', unc = FALSE, use_bam = TRUE)
        print(i, taxonIdx)
        list(po$pred_occ, ppa$pred_pot, ppl$pred_pot)
    }

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

critArith <- critLogArith <- array(0, c(length(taxa_to_fit), length(k_occ_cv), length(k_pot_cv)))
dimnames(critArith)[[1]] <- dimnames(critLogArith)[[1]] <- taxa_to_fit
dimnames(critArith)[[2]] <- dimnames(critLogArith)[[2]] <- k_occ_cv
dimnames(critArith)[[3]] <- dimnames(critLogArith)[[3]] <- k_pot_cv

for(taxonIdx in seq_along(taxa_to_fit)) {
    ## extract raw data (again) for the taxon
    taxon <- taxa_to_fit[taxonIdx]
    sub <- cell_full %>% filter(level3s == taxon)
        
    y <- sub$avg*sub$points_occ/sub$points_total  ## actual average biomass over all cells

    critArith[taxonIdx, , ] <- calc_cv_criterion(pred_occ[taxonIdx, , ], pred_pot_arith[taxonIdx, , ],
                                                 sub$points_total, y, cv_max_biomass)
    critLogArith[taxonIdx, , ] <- calc_cv_criterion(pred_occ[taxonIdx, , ], pred_pot_larith[taxonIdx, , ],
                                                    sub$points_total, y, cv_max_biomass)
}


save(critArith, critLogArith, pred_occ, pred_pot_arith, pred_pot_larith, file = file.path(interim_results_dir, 'cv_taxon_biomass.Rda'))


if(use_mpi) closeCluster(cl)
