## Extract tree-level data for all FIA plots satisfying the PalEON criteria,
## including that the data is from the most recent of the surveys for a given plot

library(dplyr)
library(readr)

## see here for FIA metadata
## https://www.fia.fs.fed.us/library/database-documentation/index.php

## Enforce column type information. See notes in fia_cols.R for why this is important.
source(file.path('preprocessing','fia_cols.R'))

## Variables to keep for our analyses:
vars = c('PLT_CN', 'PLOT', 'UNITCD', 'COUNTYCD', 'STATECD', 'SUBP', 'LAT', 'LON', 'INVYR', 
  'TREE', 'SPCD', 'DIA', 'STATUSCD', 'TPA_UNADJ', 'DRYBIO_AG')
fia = data.frame(matrix(NA, nrow = 0, ncol = length(vars)))

for(state in states){
    cat("Processing state ", state, " ... ")

    state_plot <- read_csv(file.path(raw_data_dir, paste0(state, "_PLOT.csv")),
                           col_types = cols_plot, progress = FALSE)

    state_cond <- read_csv(file.path(raw_data_dir, paste0(state, "_COND.csv")),
                           col_types = cols_cond, progress = FALSE)
    state_cond <- state_cond %>% dplyr::select(PLT_CN, STDORGCD, CONDID, COND_STATUS_CD, CONDPROP_UNADJ)

    if(exclude_plantation) {
        homogeneous_forested_plots <- state_plot %>% left_join(state_cond, by = c('CN' = 'PLT_CN')) %>%
            group_by(CN) %>% summarize(n = n(), avg_status = mean(COND_STATUS_CD), avg_origin = mean(STDORGCD)) %>%
            filter(avg_status == 1, avg_origin == 0)
    } else homogeneous_forested_plots <- state_plot %>% left_join(state_cond, by = c('CN' = 'PLT_CN')) %>%
               group_by(CN) %>% summarize(n = n(), avg_status = mean(COND_STATUS_CD)) %>%
               filter(avg_status == 1)
    
    ## only sampled, forested plots with surveys since changeover in survey design
    state_plot <- state_plot %>% filter(PLOT_STATUS_CD == 1 & INVYR >= earliest_fia_year & INVYR <= latest_fia_year)

    ## Only homogeneous condition plots; most mixed plots seem to have
    ## various amounts of non-forested

    ## Note: if we want to use plots that are partially forested, we need to retain 'COND:::CONDPROP_UNADJ'
    ## and then divide biomass/density at plot level by that number for correct scaling to area surveyed.
    ## see Github issue #2

    state_plot <- state_plot %>% filter(CN %in% homogeneous_forested_plots$CN)
    
    ## PLOTxSTATECDxUNITCDxCOUNTYCD should be a unique plot identifier
    ## select plot info only for most recent survey
    state_plot <- state_plot %>% group_by(PLOT, UNITCD, COUNTYCD) %>%
        mutate(INVYR_MAX = max(INVYR)) %>%
        filter(INVYR == INVYR_MAX)
    
    state_tree <- read_csv(file.path(raw_data_dir, paste0(state, "_TREE.csv")),
                           col_types = cols_tree, progress = FALSE)
    
    ## Only live trees with non-NA diameters, then match to plot info
    ## note that only CN is needed for matching but use of others prevents duplicated fields in result
    state_recent <- state_tree %>% filter(STATUSCD == 1 & !is.na(DIA)) %>%
        inner_join(state_plot, by = c("PLT_CN" = "CN", "PLOT"="PLOT", "INVYR" = "INVYR",
                                      "UNITCD" = "UNITCD", "COUNTYCD" = "COUNTYCD",
                                      "STATECD" = "STATECD"))
    
    state_recent <- state_recent[ , vars]
    
    cat(nrow(state_recent), " live trees in the interval ", earliest_fia_year, "-", latest_fia_year, ".\n")
    fia <- rbind(fia, state_recent)
}

assert_that(nrow(fia) == 1097212, msg = "unexpected number of trees before exclusions")

nNA <- sum(is.na(fia$SPCD))
if(nNA)
    warning("Found ", nNA, " missing taxon codes.")

## Only trees > 8 inches DBH
fia <- fia %>% filter(DIA >= diameter_cutoff_inches)
fia <- fia %>% mutate(DIA_CM = DIA / cm_to_inch)

## add PalEON level3s taxon info
fia_to_l3a <- read_csv(file.path(conversions_data_dir, fia_to_level3a_file)) %>%
    dplyr::select('fia_spcd','level3a')
l3a_to_l3s <- read_csv(file.path(conversions_data_dir, level3a_to_level3s_file)) %>%
    dplyr::select('level3a','level3s')
fia <- fia %>% left_join(fia_to_l3a, by = c("SPCD" = "fia_spcd")) %>%
    left_join(l3a_to_l3s, by = c("level3a" = "level3a"))

## Too few of certain taxa to treat separately
fia <- fia %>% mutate(level3s = ifelse(level3s %in% excluded_level3s, "Other hardwood", level3s))
           
assert_that(nrow(fia) == 399316, msg = "unexpected number of trees after exclusions")

save(fia, file = file.path(interim_results_dir, 'full_trees.Rda'))

