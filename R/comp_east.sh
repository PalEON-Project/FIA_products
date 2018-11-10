#!/bin/bash
sed -i "s/composition_region=\"western\"/composition_region=\"eastern\"/" config
Rscript -e "source('master.R'); source('stat_modeling/3_fit_composition.R')" 2>&1 3_fit_comp_east.Rout
