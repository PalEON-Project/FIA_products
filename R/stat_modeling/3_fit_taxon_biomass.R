## Fit statistical model to smooth the raw plot biomass.
## The model fits in two parts - first the proportion of points occupied by trees
## (this is much more important for the taxon-level fitting)
## then the average biomass for occupied points (called potential biomass).
## Estimated biomass is the product of occupancy and potential.

## Run-time (without cross-validation) approximately ??? with k value of 1000.

if(!exists('k_occ_taxon'))
    stop("Must specify 'k_occ_taxon'")
if(!exists('k_pot_taxon'))
    stop("Must specify 'k_pot_taxon'")

l3a_to_l3s <- read_csv(file.path(conversions_data_dir, level3a_to_level3s_file)) 
taxa_to_fit <- unique(l3a_to_l3s$level3s)
taxa_to_fit <- taxa_to_fit[!taxa_to_fit %in% excluded_level3s_OH]
print(taxa_to_fit)

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

if(do_cv) {
    stop("cv not yet set up")

} else {
    biomass_taxon <- list()
    length(biomass_taxon) <- length(taxa_to_fit)
    names(biomass_taxon) <- taxa_to_fit
    ## fit stats model
    for(taxonIdx in seq_along(taxa_to_fit)) {
        taxon <- taxa_to_fit[taxonIdx]
        sub <- cell_full %>% filter(level3s == taxon)
        k_occ <- k_occ_taxon
        k_pot <- k_pot_taxon
        ncells <- sum(sub$avg > 0)
        if(ncells < k_pot + 200)
            k_pot <- round(ncells*0.9)
        biomass_taxon[[taxonIdx]] <- try(fit(sub, newdata = pred_grid_west, k_occ = k_occ, k_pot = k_pot, return_model = TRUE, unc = TRUE, type_pot = 'log_arith', num_draws = n_stat_samples, save_draws = TRUE))
        print(taxon)
        save(biomass_taxon, file = file.path(interim_results_dir, 'fitted_taxon_biomass.Rda'))
    }
}


