## Estimate tree-level biomass using diameter of trees and PEcAn allometries

## Run-time: approximately 3 hours

library(dplyr)
library(readr)
library(PEcAn.allometry)
library(PEcAn.logger)

load(file.path(interim_results_dir, 'full_trees.Rda'))

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

## This lists components: we are using 6 (stem) because it has the most data, though 
## would like to use component 2 (Whole tree (aboveground)) or component 3
## (Whole tree (above stump))
if(FALSE) {
    data(allom.components)
    allom.components
}

## Fit allometry models for all taxa
unique_pecan_allom <- unique(fia$pecan_allometry_spcd)
pecan_taxa <- lapply(seq_len(length(unique_pecan_allom)), function(i) 
    data.frame(spcd = as.numeric(strsplit(unique_pecan_allom[i], split = ',')[[1]])))
names(pecan_taxa) <- unique_pecan_allom

## first check which allometries are already in the allom_dir and don't redo
## TBD
allom_stats = load.allom(allom_dir)

## create via loop rather than passing list of dataframes because of PEcAn bug -
## errors out for spcd=701
for(i in seq_along(pecan_taxa)) {
    if(!names(pecan_taxa)[i] %in% names(allom_stats)) {
        cat("Fitting ", names(pecan_taxa)[i], ".\n")
        allom_stats[[names(pecan_taxa)[i]]] <- try(AllomAve(pecan_taxa[i], ngibbs=1000, components = 6,
                             outdir = allom_dir, dmin = 10, dmax = 150))
}}

## tricky cases from initial use of fia_to_pecan_v0.3.csv:
## one allometry, lmtd dbh range: 763, 920, 931
## one allometry: 126, 332, 355, 543, 544, 691, 743, 762, 970
## two allometry: 93, 241, 731, 827, 835, 901
## no fitting done: 10, 68, 401, 407, 491, 540, 711, 741, 809, 823, 972 (no component 6)
## fit failures with R error: 315, 319, 701 - these are because max diam less than cutoff
## 931 (sassafras) is highly uncertain and giving crazy biomasses, presumably because of limited data and dbh range in (11,14)
## 355 (serviceberry) gives some outliers
## 762 (black cherry) gives some outliers
## 371 (yellow birch) has some outliers but fair amount of allometries

## not_fit <- c(10, 68, 401, 407, 491, 540, 711, 741, 809, 823, 972)
## fit_failure <- c(315, 319, 701)
## fit_unstable <- c(355, 762, 763, 931)
## missing_pecan_taxa <- c(not_fit, fit_failure, fit_unstable)

## tricky cases from fia_to_pecan_v0.4.csv (first pass at v0.4):
## one allometry: 970;972;974;975;977 (2879 trees),  760;762 (9280 trees),
## two allometries: 315;471;491;319;355;356;357;931;935;391;701 (657 trees),  315;471;491;319;355;356;357;931;935 (2065 trees)
## no fitting done: 68;66;67 (827 trees), 600;601;602 (1611 trees)

## These cases from first pass at 0.4 seem to have been dealt with in final v0.4.

## site effects seem to product crazy allometries - mu0 and mu1 have outliers and taus can be big; how many allometries do we need for stability?

ntrees <- nrow(fia)

tmp_plt <- rep("", ntrees)
tmp_subplt <- rep(0, ntrees)
tmp_tree <- rep(0, ntrees)
biomass <- matrix(0, nrow = ntrees, ncol = n_allom_samples)

counter <- 1
## We want to run the allometry prediction for all trees of a given taxon in a given plot simultaneously because this allows for shared allometric parameters (but different tree effects) for trees of the same taxon. In this case we simply run the prediction for all trees in a plot at once.
## TREE is unique only within subplot.
## This takes 3.5 hours on smeagol on one core.
set.seed(1)
for(plt in unique(fia$PLT_CN)) {
    sub <- fia %>% filter(PLT_CN == plt) %>% dplyr::select(PLT_CN, SUBP, TREE, DIA_CM, pecan_allometry_spcd)
    pred <- allom.predict(allom_dir,
                          dbh = sub$DIA_CM,
                          pft = sub$pecan_allometry_spcd,
                          component = allom_component, 
                          use = 'Bg', # 'mu' in unstable statistically
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

save(tmp_plt, tmp_tree, tmp_subplt, biomass,
     file = file.path(interim_results_dir('biomass_samples.Rda'))

## need to check back as to whether we want to do this
if(FALSE) 
    biomass[biomass > biomass_max_kg] <- biomass_max_kg

if(!do_allom_uncertainty) {
    tmp = rowMeans(biomass)
    fia <- fia %>% mutate(biomass = tmp)
} else {
    biomass_tbl <- cbind(data_frame(PLT_CN = tmp_plt, SUBP = tmp_subplt, TREE = tmp_tree), biomass)
    names(biomass_tbl)[4:ncol(biomass_tbl)] <- paste0('biomass_smp_', 1:n_allom_samples)

    fia <- fia %>% inner_join(biomass_tbl, by = c("PLT_CN", "SUBP", "TREE"))
    if(nrow(fia) != ntrees)
        stop("Something wrong with merge of biomass samples with tree data")
}

save(fia, file = file.path(interim_results_dir, 'full_trees_with_biomass.Rda'))


