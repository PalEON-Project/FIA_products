#!/bin/bash
sed -i "s/composition_region=\"western\"/composition_region=\"eastern\"/" config
sed -i "s/first150k/second150k/" config
Rscript -e "source('master.R'); source('composition_modeling/3_draw_proportions.R')" 2>&1 draw_props_eastern_second150k.Rout
