#!/bin/bash
cd /storage2/asirix/Degron_aux
wget wget -P /. http://172.26.5.30/ftp/by_Project/Degron_3D/Degrone_Hi-C/2022-06-18-degron-main-from-Kostya/cohesin_condensin_hic_only/RA1_30.hic
wget wget -P /. http://172.26.5.30/ftp/by_Project/Degron_3D/Degrone_Hi-C/2022-06-18-degron-main-from-Kostya/cohesin_condensin_hic_only/RC1_30.hic
for y in 5000 25000
do
hic2cool convert -r $y -p 10 ./RA1_30.hic ./'RA1_'$y'.cool'
hic2cool convert -r $y -p 10 ./RC1_30.hic ./'RC1_'$y'.cool'
cooler balance ./'RA1_'$y'.cool'
cooler balance ./'RC1_'$y'.cool'
done 
for i in 5 8 9
do
coolpup.py ./RC1_5000.cool ./'TADs_min_ins_0.'$i'.bed' --rescale --local --nproc 0
coolpup.py ./RC1_25000.cool ./'TADs_min_ins_0.'$i'.bed' --rescale --local --nproc 0
coolpup.py ./RA1_5000.cool ./'TADs_min_ins_0.'$i'.bed' --rescale --local --nproc 0
coolpup.py ./RA1_25000.cool ./'TADs_min_ins_0.'$i'.bed' --rescale --local --nproc 0 
done
for g in 5 8 9 
do
plotpup.py --dpi 600 --input_pups ./'RA1_5000.cool-5.0K_over_TADs_min_ins_0.'$g'_10-shifts_local_rescaled.clpy' --no_score --not_symmetric --output ./'RA1_5k_0.'$g'.png'
plotpup.py --dpi 600 --input_pups ./'RA1_25000.cool-25.0K_over_TADs_min_ins_0.'$g'_10-shifts_local_rescaled.clpy' --no_score --not_symmetric --output ./'RA1_25k_0.'$g'.png'
plotpup.py --dpi 600 --input_pups ./'RC1_5000.cool-5.0K_over_TADs_min_ins_0.'$g'_10-shifts_local_rescaled.clpy' --no_score --not_symmetric --output ./'RC_5k_0.'$g'.png'
plotpup.py --dpi 600 --input_pups ./'RC1_25000.cool-25.0K_over_TADs_min_ins_0.'$g'_10-shifts_local_rescaled.clpy' --no_score --not_symmetric --output ./'RC1_25k_0.'$g'.png'
done
