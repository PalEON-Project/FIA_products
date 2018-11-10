#!/bin/bash
sed -i "s/composition_region=\"eastern\"/composition_region=\"western\"/" config
sed -i "s/first150k/second150k/" config
Rscript -e "source('master.R'); source('composition_modeling/3_draw_proportions.R')" 2>&1 draw_props_western_second150k.Rout
