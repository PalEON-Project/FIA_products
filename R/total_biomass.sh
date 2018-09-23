#!/bin/bash
Rscript -e "source('master.R'); source('biomass_modeling/2_fit_total_biomass.R')" 2>&1 2_fit_total_biomass.Rout
