#!/bin/bash
Rscript -e "source('master.R'); source('biomass_modeling/2_fit_taxon_biomass.R')" 2>&1 2_fit_taxon_biomass.Rout
