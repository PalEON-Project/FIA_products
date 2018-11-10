#!/bin/bash
Rscript -e "source('master.R'); source('density_modeling/1_cv_total_density.R')" 2>&1 1_cv_total_density.Rout
