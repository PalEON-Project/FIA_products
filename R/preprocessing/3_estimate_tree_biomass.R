## Estimate tree-level biomass using diameter of trees and PEcAn allometries

## Run-time: This takes 3.5 hours for 10 samples each on modern desktop on one core.
## Most computational time is spent in predicting biomass for the individual FIA trees.

library(dplyr)
library(assertthat)
library(readr)
library(PEcAn.allometry)
library(PEcAn.logger)

load(file.path(interim_results_dir, 'full_trees.Rda'))

## Read in conversion info to go from FIA taxa to PalEON-determined aggregations of taxa
## for which we can estimate statistically stable allometric relationships.


if(use_agb) {
    cols_conv <- cols(
        level3a = col_character(),
        pecan_allometry_spcd = col_character(),  ## multiple spcd separated by ; for non-exact matches
        pecan_allometry_common_names = col_character()
    )

    ## Note: will see warning about missing column names for X4-X12
    taxa_conversion <- read_csv(file.path(conversions_data_dir, paleon_to_chojnacky_conversion_file),
                                col_types = cols_conv) %>%
        dplyr::select(level3a, beta0, beta1)
    

    fia <- fia %>% left_join(taxa_conversion, by = c("level3a" = "level3a")) %>%
        rename(int = beta0, slope = beta1)

    fia$int[fia$level3a == 'Unknown tree'] <- taxa_conversion$beta0[taxa_conversion$level3a == "Other hardwood"]
    fia$slope[fia$level3a == 'Unknown tree'] <- taxa_conversion$beta1[taxa_conversion$level3a == "Other hardwood"]
    fia$int[fia$level3a == 'Douglas fir'] <- taxa_conversion$beta0[taxa_conversion$level3a == "Fir"]
    fia$slope[fia$level3a == 'Douglas fir'] <- taxa_conversion$beta1[taxa_conversion$level3a == "Fir"]
    
    predict_biomass <- function(dbh_cm, b0, b1) {
        return(exp(b0+b1*log(dbh_cm)))
    }

    fia <- fia %>% mutate(biomass = predict_biomass(DIA_CM, int, slope)) 

    assert_that(min(fia$biomass) > 0 & max(fia$biomass) < 1e6, msg = "extreme biomass values")

} else {
    cols_conv <- cols(
        fia_spcd = col_integer(),
        fia_counts = col_integer(),
        fia_scientific_names = col_character(),
        fia_common_names = col_character(),
        pecan_allometry_spcd = col_character(),  ## multiple spcd separated by ; for non-exact matches
        pecan_allometry_common_names = col_character(),
        acronym = col_character(),
        pft = col_character()
    )
    taxa_conversion <- read_csv(file.path(conversions_data_dir, fia_to_pecan_conversion_file),
                            col_types = cols_conv) %>%
    dplyr::select(fia_spcd, pecan_allometry_spcd)
                               

    ## only taxa with a PEcAn allometry (hopefully this will include all trees in dataset)
    taxa_conversion <- taxa_conversion %>% filter(!is.na(pecan_allometry_spcd)) %>%
        mutate(pecan_allometry_spcd = gsub(";", ",", pecan_allometry_spcd))
    
    
    fia <- fia %>% left_join(taxa_conversion, by = c("SPCD" = "fia_spcd")) 
    
    assert_that(sum(is.na(fia$pecan_allometry_spcd)) == 0, msg = "missing allometries")
    
    ## Fit allometry models for all taxa
    
    unique_pecan_allom <- unique(fia$pecan_allometry_spcd)
    pecan_taxa <- lapply(seq_len(length(unique_pecan_allom)), function(i) 
        data.frame(spcd = as.numeric(strsplit(unique_pecan_allom[i], split = ',')[[1]])))
    names(pecan_taxa) <- unique_pecan_allom
    
    ## First check which allometries are already in the allom_dir and don't redo
    allom_stats = load.allom(allom_dir)
    
    ## Fit allometry, creating via loop rather than passing list of dataframes because of PEcAn bug
    for(i in seq_along(pecan_taxa)) {
        if(!names(pecan_taxa)[i] %in% names(allom_stats)) {
            cat("Fitting ", names(pecan_taxa)[i], ".\n")
            allom_stats[[names(pecan_taxa)[i]]] <- try(AllomAve(pecan_taxa[i], ngibbs=1000,
                                                                components = allom_component,
                                                                outdir = allom_dir, dmin = dmin_allom_fit, dmax = dmax_allom_fit))
        }}
    
    ## Setup storage.
    n_trees <- nrow(fia)
    tmp_plt <- rep("", n_trees)
    tmp_subplt <- rep(0, n_trees)
    tmp_tree <- rep(0, n_trees)
    biomass <- matrix(0, nrow = n_trees, ncol = n_allom_samples)
    
    
    ## We want to run the allometry prediction for all trees of a given taxon in a given plot simultaneously
    ## because this allows for shared allometric parameters (but different tree effects) for trees of the
    ## same taxon. In this case we simply run the prediction for all trees in a plot at once.
    ## TREE is unique only within subplot.
    set.seed(1)
    counter <- 1
    for(plt in unique(fia$PLT_CN)) {
        sub <- fia %>% filter(PLT_CN == plt) %>% dplyr::select(PLT_CN, SUBP, TREE, DIA_CM, pecan_allometry_spcd)
        ## Inclusion of site effects seem to product crazy allometries:
        ## mu0 and mu1 have outliers and taus can be big, so have 'use' be 'Bg'
        pred <- allom.predict(allom_dir,
                              dbh = sub$DIA_CM,
                              pft = sub$pecan_allometry_spcd,
                              component = allom_component, 
                              use = 'Bg', # 'mu' is unstable statistically
                              n = n_allom_samples,
                              interval = "prediction", 
                              single.tree = FALSE)
        nsub <- nrow(sub)
        tmp_plt[counter:(counter+nsub-1)] <- sub$PLT_CN
        tmp_subplt[counter:(counter+nsub-1)] <- sub$SUBP
        tmp_tree[counter:(counter+nsub-1)] <- sub$TREE
        biomass[counter:(counter+nsub-1), ] <- t(pred)
        counter <- counter + nsub
        print(plt)
    }

    assert_that(min(biomass) > 0 & max(biomass) < 1e6, msg = "extreme biomass values")

    biomass[biomass > biomass_max_kg] <- biomass_max_kg
    
    save(tmp_plt, tmp_tree, tmp_subplt, biomass,
         file = file.path(interim_results_dir, 'biomass_samples.Rda'))
} 
 

if(!do_allom_uncertainty) {
    fia <- fia %>% mutate(biomass = rowMeans(biomass))
} else {
    biomass_tbl <- cbind(data_frame(PLT_CN = tmp_plt, SUBP = tmp_subplt, TREE = tmp_tree), biomass)
    names(biomass_tbl)[4:ncol(biomass_tbl)] <- paste0('biomass_smp_', 1:n_allom_samples)
    
    fia <- fia %>% inner_join(biomass_tbl, by = c("PLT_CN", "SUBP", "TREE"))
    assert_that(nrow(fia) == n_trees,
                msg = "Something wrong with merge of biomass samples with tree data")
}

fn <- 'full_trees_with_biomass'
if(use_agb)
    fn <- paste0(fn, '_agb')
fn <- paste0(fn, '.Rda')

save(fia, file = file.path(interim_results_dir, fn))


    

