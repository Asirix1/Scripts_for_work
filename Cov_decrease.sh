#!/bin/bash
for i in 01 05 1 5 1 2 3 5 6:
do
cd /storage2/asirix/Degron_aux/Cov_decrease/RA1/
mkdir ./'Reads_ratio_'$i
cd ./'Reads_ratio_'$i
awk 'BEGIN  {srand()} !/^$/  { if (rand() <= .$i || FNR==1) print > "merged_nodups_sample.txt"}' $1
awk -v OFS="\t" '{$1=$1; print}' ./merged_nodups_sample.txt > ./merged_nodups_sample_tab.txt   
cooler cload pairs -c1 2 -p1 3 -c2 6 -p2 7 /storage2/asirix/Degron_aux/Cov_decrease/mm10.chr.sizes:5000 ./merged_nodups_sample_tab.txt ./Decrease.cool
cooler balance ./Decrease.cool
coolpup.py ./Decrease.cool /storage2/asirix/Degron_aux/TADs_min_ins_0.0.bed --rescale --rescale_flank 50 --local --nproc 0 -o ./Cov_decr.plpy
plotpup.py --dpi 600 --input_pups ./Cov_decr.plpy --output ./Cov_decr.png
done