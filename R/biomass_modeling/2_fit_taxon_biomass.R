## Fit statistical model to smooth the raw plot biomass.
## The model fits in two parts - first the proportion of points occupied by trees
## (this is much more important for the taxon-level fitting)
## then the average biomass for occupied points (called potential biomass).
## Estimated biomass is the product of occupancy and potential.

## Run-time: this takes about 3 hours, if done in parallel across the 21 taxa.

load(file.path(interim_results_dir, 'full_trees_with_biomass_grid.Rda'))

if(!exists('k_occ_taxon'))
    stop("Must specify 'k_occ_taxon'")
if(!exists('k_pot_taxon'))
    stop("Must specify 'k_pot_taxon'")

taxa_to_fit <- unique(fia$level3s)
taxa_to_fit <- taxa_to_fit[!is.na(taxa_to_fit)]  ## 22 Douglas fir trees in dataset. Not fit.
print(taxa_to_fit)

## Allow for parallelization across taxa, including on Berkeley Statistics cluster with SLURM scheduler
library(doParallel)
if(n_cores == 0) {
    if(Sys.getenv("SLURM_JOB_ID") != "") {
        n_cores <- Sys.getenv("SLURM_CPUS_PER_TASK")
    } else n_cores <- detectCores()
}

registerDoParallel(cores = n_cores)

## for naming consistency with PLS, will refer to FIA plots as 'points'

## total number of fia points per grid cell
plots_per_cell <- fia %>% group_by(PLT_CN, x, y, cell) %>%
    summarize(n_trees = n()) %>% 
    group_by(cell) %>% summarize(points_total = n()) 

area_conversion <- acres_per_ha * squ_feet_per_acre /
    (pi * fia_subplot_radius_in_feet^2 * fia_subplots_per_plot) / kg_per_Mg 

## Total biomass per taxon per plot only for taxa represented in a plot.
## Currently grouping taxa by level3s but this could change.
biomass_taxon_plot <- fia %>% group_by(PLT_CN, level3s, cell, x, y) %>%
    summarize(total_biomass = sum(biomass) * area_conversion)
## after conversion, biomass is Mg/Ha

## Average biomass per taxon per cell only for taxa represented in a cell,
## and only for plots in which taxon is present, so that potential biomass
## modeling ignores occupancy.
biomass_taxon_cell <- biomass_taxon_plot %>% group_by(cell, x, y, level3s) %>%
    summarize(avg = mean(total_biomass),
              geom_avg = mean(log(total_biomass)), points_occ = n()) 

## Flesh out data so have 0s for taxa missing in each cell
## (needed for occupancy part of model)
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

assert_that(nrow(cell_full) == nrow(plots_per_cell) * length(taxa_to_fit),
            msg = "Dataset does not have data for all taxa for all occupied cells.")

## Fit stats model to each taxon, in parallel
biomass_taxon <- foreach(taxonIdx = seq_along(taxa_to_fit)) %dopar% {
    taxon <- taxa_to_fit[taxonIdx]
    sub <- cell_full %>% filter(level3s == taxon)
    ## try() allows fitting to fail for a taxon without stopping fitting for other taxa
    try(fit(sub, newdata = pred_grid_paleon, k_occ = k_occ_taxon, k_pot = k_pot_taxon,
            return_model = TRUE, unc = TRUE, type_pot = 'arith', num_draws = n_stat_samples,
            save_draws = TRUE, use_bam = TRUE))
}

names(biomass_taxon) <- taxa_to_fit
save(biomass_taxon, file = file.path(interim_results_dir, 'fitted_taxon_biomass.Rda'))
    



