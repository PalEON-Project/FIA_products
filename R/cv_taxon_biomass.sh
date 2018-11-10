#!/bin/bash
Rscript -e "source('master.R'); source('biomass_modeling/1_cv_taxon_biomass.R')" 2>&1 1_cv_taxon_biomass.Rout
