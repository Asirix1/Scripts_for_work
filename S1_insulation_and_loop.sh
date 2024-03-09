#!/bin/bash
for i in 
do
conda activate aapopov
cd /storage3/asirix/S1/Insulation/merged_nodups/
gunzip 'storage3/asirix/S1/Insulation/merged_nodups/'$i
i=${i::-3}
cooler cload pairs -c1 2 -p1 3 -c2 6 -p2 7 /storage3/asirix/S1/Insulation/reference.chromsizes.for_cooler:5000 'storage3/asirix/S1/Insulation/merged_nodups/'$i '/storage3/asirix/S1/Insulation/Cool/'$i'.cool'
cooler balance '/storage3/asirix/S1/Insulation/Cool/'$i'.cool'
coolpup.py '/storage3/asirix/S1/Insulation/Cool/'$i'.cool' /storage3/asirix/S1/Insulation/Tracks/ctcf_pairs.bedpe --flank 200000 --nproc 0 -o '/storage3/asirix/S1/Insulation/Loop/'$i'_loop.clpy'
plotpup.py --dpi 600 --input_pups '/storage3/asirix/S1/Insulation/Loop/'$i'_loop.clpy' --output '/storage3/asirix/S1/Insulation/Figures/'$i'_loop.png'
coolpup.py '/storage3/asirix/S1/Insulation/Cool/'$i'.cool' /storage3/asirix/S1/Insulation/Tracks/ENCFF660GHM_peaks.bed --local --flank 200000 --nproc 0 -o '/storage3/asirix/S1/Insulation/Insulation/'$i'_insulation.clpy'
plotpup.py --dpi 600 --input_pups '/storage3/asirix/S1/Insulation/Insulation/'$i'_insulation.clpy' --output '/storage3/asirix/S1/Insulation/Insulation/'$i'_insulation.png'
done