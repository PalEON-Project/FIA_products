## Full workflow for FIA analyses
## Chris Paciorek
## 2018-02-01

## note required input files (download from PalEON Wiki)

# 7z -e NRS_actual_coords_forested_plots.zip
# unzip NRS_actual_coords_forested_plots_2nd_run.zip
# unzip NRS_actual_coords_forested_plots_3rd_run.zip


## setup R packages (with version control) 

if(FALSE) {  ## turn off during initial development  
    checkpoint::checkpoint("2018-01-18")
    library(devtools)
    install_github("PecanProject/pecan",subdir="base/logger", ref = pecan_base_logger_commit)
    install_github("PecanProject/pecan",subdir="modules/allometry", ref = pecan_modules_allometry_commit)
}

## get configuration variables

source('config')



## setup directories and basic data

states = c('MN', 'WI', 'MI', 'IL', 'IN', 'OH', 'PA', 'NY', 'NJ', 'ME', 'VT', 'MA', 'CT', 'NH', 'RI')

paleon_regions_west <- c(2,3,5,6,11,12) ## MI/IN and west
paleon_regions_west_ohio <- c(2,3,5,6,9,11,12) ## OH and west
paleon_regions_east <- c(1,4,7,8,9,10)
paleon_regions <- 1:12

paleon_states_west <- c(17, 18, 26, 27, 55) 
paleon_states_west_ohio <- c(paleon_states_west, 39)
paleon_states_east <- c(9, 23, 25, 33, 34, 36, 39, 42, 44, 50)  # includes Ohio (39)


## setup directories relative to current working directory

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

## source files with R functions
code_files <- list.files(code_dir, pattern = ".R$", full.names = TRUE)
sapply(code_files, source)

excluded_level3s_OH <- c('Dogwood', 'Chestnut', 'None')
## Dogwood: only 20 trees >= 20 cm
## Chestnut: no modern data

## values used in cross-validation for biomass/density
k_occ_cv <- c(100,250,500,1000,1500,2000,2500)  # 3000,3500
k_pot_cv = c(100,250,500,1000,1500,2000,2500,3000,3500)

## TMP: for seeing all cols instead of tibble
adf <- as.data.frame
