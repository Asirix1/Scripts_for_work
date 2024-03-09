#!/bin/bash
source /mnt/storage/home/aapopov/anaconda3/bin/activate
conda activate /mnt/storage/home/aapopov/anaconda3/envs/aapopov
for i in DNase.bed ATAC.bed
do
export i
for g in Offset1_MNase304U.coverage.bw Offset1_MNase79.9U.coverage.bw Offset1_MNase20.6U.coverage.bw Offset1_MNase5.4U.coverage.bw Offset1_MNase304U.coverage.bw Offset1_MNase79.9U.coverage.bw Offset1_MNase20.6U.coverage.bw Offset1_MNase5.4U.coverage.bw newS1-Offset1-bin1.coverage.bw K562S1-Offset1-bin1.bigwig K562MNase-Offset1-bin1.bigwig K562DNase-Offset1-bin1.bigwig Offset1_MNase304U.coverage.bw Offset1_MNase79.9U.coverage.bw Offset1_MNase20.6U.coverage.bw Offset1_MNase5.4U.coverage.bw newS1-Offset1-bin1.coverage.bw K562WG-Offset1-bin1.bigwig K562S1-Offset1-bin1.bigwig K562MNase-Offset1-bin1.bigwig K562DNase-Offset1-bin1.bigwig 
do
qsub <<- BODY
#PBS -l select=1:ncpus=24:mem=192g,walltime=1:0:0 -q dl560g10q@vm-pbs2
#PBS -N aapopov 
export g
python /mnt/storage/home/aapopov/S1_final/Drawing_graphs_S1.py
BODY
done 
done
