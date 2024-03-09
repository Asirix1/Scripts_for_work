#!/bin/bash
for y in 10000 25000 100000
do
cd /storage2/asirix/Degron_aux
mkdir ./'Resolution_'$y
cd ./'Resolution_'$y
hic2cool convert -r $y -p 10 /storage2/asirix/Degron_aux/RA1_30.hic ./'RA1_'$y'.cool'
hic2cool convert -r $y -p 10 /storage2/asirix/Degron_aux/RC1_30.hic ./'RC1_'$y'.cool'
cooler balance ./'RA1_'$y'.cool'
cooler balance ./'RC1_'$y'.cool'
m=50
for i in 0 5 8 9
do
coolpup.py /storage2/asirix/Degron_aux/'Resolution_'$y/'RA1_'$y'.cool' /storage2/asirix/Degron_aux/'TADs_min_ins_0.'$i'.bed' --rescale --rescale_flank $m --local --nproc 0
coolpup.py /storage2/asirix/Degron_aux/'Resolution_'$y/'RC1_'$y'.cool' /storage2/asirix/Degron_aux/'TADs_min_ins_0.'$i'.bed' --rescale --rescale_flank $m --local --nproc 0
plotpup.py --dpi 600 --input_pups ./RA1_*.clpy --no_score --output ./'RA1_0.'$i'.png'
plotpup.py --dpi 600 --input_pups ./RC1_*.clpy --no_score --output ./'RC1_0.'$i'.png'
done
done
