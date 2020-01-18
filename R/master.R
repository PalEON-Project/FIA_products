## Full workflow for FIA analyses
## Chris Paciorek
## 2019-06-02 (original 2018-02-01)
## Analyses done using R 3.5.1
## as shown below, uses various R packages from 2019-04-13


## Various required input files (download from PalEON Wiki)

# 7z -e NRS_actual_coords_forested_plots.zip
# unzip NRS_actual_coords_forested_plots_2nd_run.zip
# unzip NRS_actual_coords_forested_plots_3rd_run.zip

## Taxon conversion files

## Based on discussions with Charlie Cogbill, Mike Dietze:
## Convert FIA taxa to PEcAn allometry taxa: fia_to_pecan_v0.4.csv
## when using stem biomass (allometry component 6)

## But as of May/June 2019 we are focusing on the Chojnacky
## AGB allometry
##
## level3s is the aggregation for most of our analyses,
## though we do exclude some uncommon taxa.

## get configuration variables
source('config')


## setup R packages (with version control) 
if(use_original_R_package_versions) {  
    ## bug in checkpoint: errors out with 'arguments imply differing number of rows':
    ## work-around is to avoid wget as download method
    if(options('download.file.method') == 'wget')
        options('download.file.method' = NULL)
    checkpoint::checkpoint(checkpoint_date, R.version = R_version)
    devtools::install_github("PecanProject/pecan",subdir="base/logger", ref = pecan_base_logger_commit)
    devtools::install_github("PecanProject/pecan",subdir="modules/allometry", ref = pecan_modules_allometry_commit)
}


## geographic subsetting

states = c('MN', 'WI', 'MI', 'IL', 'IN', 'OH', 'PA', 'NY', 'NJ', 'ME', 'VT', 'MA', 'CT', 'NH', 'RI')

## PalEON-defined regions: some PalEON input files are stratified
## by our own regions that cross state boundaries
paleon_regions_west <- c(2,3,5,6,11,12) ## MI/IN and west
paleon_regions_west_ohio <- c(2,3,5,6,9,11,12) ## OH and west
paleon_regions_east <- c(1,4,7,8,9,10)
paleon_regions <- 1:12

## State FIPS codes:
paleon_states_west <- c(17, 18, 26, 27, 55) 
paleon_states_west_ohio <- c(paleon_states_west, 39)
paleon_states_east <- c(9, 23, 25, 33, 34, 36, 39, 42, 44, 50)  # includes Ohio (39)


## setup directories relative to current working directory,
## if directories not set in config file
if(raw_data_dir == "")
    raw_data_dir <- file.path("..", "data", "raw") 
if(conversions_data_dir == "")
    conversions_data_dir <- file.path("..", "data", "conversions") 
if(allom_dir == "")
    allom_dir <- file.path("..", "data", "allom")
if(code_dir == "")
    code_dir <- file.path("code")
if(output_dir == "")
    output_dir <- file.path("..", "output")
if(interim_results_dir == "")
    interim_results_dir <- file.path("..", "data", "interim")

## source all files with R functions
code_files <- list.files(code_dir, pattern = ".R$", full.names = TRUE)
sapply(code_files, source)

excluded_level3s <- c('Dogwood', 'Chestnut')
## Dogwood: only 3 trees >= 20 cm
## Chestnut: no modern data

## values used in cross-validation for biomass/density
## caution: k values of 3000,3500 for 'occ' can take a very long time to fit
k_occ_cv <- c(100,250,500,1000,1500,2000,2500) 
k_pot_cv <- c(100,250,500,1000,1500,2000,2500,3000,3500)


## The following documents the workflow to do the analyses
## but would not generally be run all at once,
## so embedded in conditional statements.

## One needs to source master.R before running any of these.

## for data download and processing:
if(FALSE) {
    source(file.path("preprocessing", "1_download_fia.R"))
    source(file.path("preprocessing", "2_process_for_paleon.R"))
    ## next step not needed for composition analysis
    ## but 4_setup_grid.R assumes use of output from 3_estimate_tree_biomass.R
    source(file.path("preprocessing", "3_estimate_tree_biomass.R"))
    source(file.path("preprocessing", "4_setup_grid.R"))
}

## for biomass modeling
if(FALSE) {
    ## to determine optimal k_occ, k_pot, run these two lines:
    source(file.path("biomass_modeling", "1_cv_total_biomass.R"))
    source(file.path("biomass_modeling", "1_cv_taxon_biomass.R"))
    ## based on results, set k_occ_taxon, k_pot_taxon, k_pot_total in 'config'
    ## then run:
    source(file.path("biomass_modeling", "2_fit_total_biomass.R"))
    source(file.path("biomass_modeling", "2_fit_taxon_biomass.R"))
    source(file.path("biomass_modeling", "3_output_biomass.R"))
}

## for composition modeling
if(FALSE) {
    source(file.path("composition_modeling", "1_setup_composition.R"))
    ## run next line twice, changing 'composition_region' in config
    ## to be 'eastern' and 'western'
    ## because fitting takes 6-7 days for 150k MCMC iterations
    source(file.path("composition_modeling", "2_fit_composition.R"))
    source(file.path("composition_modeling", "3_draw_proportions.R"))
    source(file.path("composition_modeling", "4_output_composition.R"))
}

## temporary: convert tibble to regular data frame for easier viewing in R
adf <- as.data.frame
