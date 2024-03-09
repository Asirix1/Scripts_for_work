#!/bin/bash
qsub <<- BODY
#PBS -l select=1:ncpus=24:mem=192g,walltime=71:0:0 -q dl560g10q@vm-pbs2
#PBS -N $1
source /opt/shared/anaconda/anaconda3-2020/bin/activate
conda activate aapopov
macs3 callpeak -t ./Galaxy19-\[Filter_SAM_or_BAM\,_output_SAM_or_BAM_on_data_17__bam\].bam --nomodel --extsize 131
macs3 callpeak -t ./Galaxy20-\[Filter_SAM_or_BAM\,_output_SAM_or_BAM_on_data_17__bam\].bam --nomodel --extsize 131
