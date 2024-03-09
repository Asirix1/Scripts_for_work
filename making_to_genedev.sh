#!/bin/bash

for y in 718 718_a Ch_a
do

cp /data/aapopov/Projects/HumanDros/Zapusk_Sirius_Human_Dros/$y'/aligned/inter.hic' /data/aapopov/Projects/HumanDros/Zapusk_Sirius_Human_Dros/Zapusk_sirius_12_23/$y/$y'.hic'
cp /data/aapopov/Projects/HumanDros/Zapusk_Sirius_Human_Dros/$y'/aligned/inter.txt' /data/aapopov/Projects/HumanDros/Zapusk_Sirius_Human_Dros/Zapusk_sirius_12_23/$y/$y'_stats.txt'
cp /data/aapopov/Projects/HumanDros/Zapusk_Sirius_Human_Dros/$y'/aligned/inter_30.hic' /data/aapopov/Projects/HumanDros/Zapusk_Sirius_Human_Dros/Zapusk_sirius_12_23/$y/$y'_30.hic'
cp /data/aapopov/Projects/HumanDros/Zapusk_Sirius_Human_Dros/$y'/aligned/inter_30.txt' /data/aapopov/Projects/HumanDros/Zapusk_Sirius_Human_Dros/Zapusk_sirius_12_23/$y/$y'_30_stats.txt'

done
