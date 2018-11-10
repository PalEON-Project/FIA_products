#!/bin/bash
sed -i "s/composition_region=\"eastern\"/composition_region=\"western\"/" config
Rscript -e "source('master.R'); source('stat_modeling/3_fit_composition.R')" 2>&1 3_fit_comp_west_run2.Rout
