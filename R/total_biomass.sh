#!/bin/bash
Rscript -e "source('master.R'); source('stat_modeling/1_setup_data.R'); source('stat_modeling/3_fit_total_biomass.R')" 2>&1 3_fit_total_biomass.Rout
