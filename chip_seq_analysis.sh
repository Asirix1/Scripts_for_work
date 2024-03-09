#!/bin/bash

for DIR in */
do
DIR=${DIR%?}
cd /data/aapopov/Projects/HumanDros/18_01_24_BGI_hi_c_CTCF/ChIP_seq_S2/$DIR
bowtie2 \
--phred33 \
--mm \
--very-sensitive \
--threads 20 \
-x /data/aapopov/Projects/HumanDros/18_01_24_BGI_hi_c_CTCF/ChIP_seq_S2/genome_drosophila/dm6 \
-1 /data/aapopov/Projects/HumanDros/18_01_24_BGI_hi_c_CTCF/ChIP_seq_S2/${DIR}/${DIR}_L1_1.fq.gz \
-2 /data/aapopov/Projects/HumanDros/18_01_24_BGI_hi_c_CTCF/ChIP_seq_S2/${DIR}/${DIR}_L1_2.fq.gz \
2> \
bowtie2.log \
| samtools view -h -b - > ${DIR}_aligned.bam

samtools view -h -b -F 4 -q 30 -@ 10 -o ${DIR}_aligned_filtered.bam ${DIR}_aligned.bam 
samtools sort -n -O BAM -@ 10 ${DIR}_aligned_filtered.bam  \
| samtools fixmate -m -@ 10 - - \
| samtools sort -O BAM -@ 10 \
| samtools markdup -r -S -@ 10 - ${DIR}_aligned_filtered_sorted_duprmv.bam
samtools index ${DIR}_aligned_filtered_sorted_duprmv.bam ${DIR}_aligned_filtered_sorted_duprmv.bam.bam.bai
samtools flagstat ${DIR}_aligned.bam > ${DIR}_aligned.bam.log
samtools flagstat ${DIR}_aligned_filtered_sorted_duprmv.bam > ${DIR}_aligned_filtered_sorted_duprmv.bam.log
 
samtools flagstat -@ 3 data/processed/CTCF/Bowtie2/CTCF_Rep1_ENCFF001HLV_trimmed_aligned_filt_sort_nodup.bam > \
analysis/AlignmentStats/ChIP/CTCF_Rep1_ENCFF001HLV_trimmed_aligned_filt_sort_nodup.flagstat.log
samtools sort -O BAM -n -@ 3 data/processed/CTCF/Bowtie2/CTCF_Rep1_ENCFF001HLV_trimmed_aligned_nofilt.bam \

samtools sort -O BAM -n -@ 3 ${DIR}_aligned.bam \
| samtools fixmate -m -@ 3 - - \
| samtools sort -O BAM -@ 3 \
| samtools markdup -@ 3 - - \
| samtools flagstat - > ${DIR}_aligned_nofilt_dupmarked_flagstat.log

bamCoverage -b ${DIR}_aligned_filtered_sorted_duprmv.bam -o ${DIR}_aligned_filtered_sorted_duprmv.bigwig
done