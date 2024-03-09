#!/bin/bash
for t in 13 14 15 16 18 19 20 22
do
gunzip ./'PICO'$t'_R1.fq.gz'
cat ./'PICO'$t'_R1.fq' > ./PICO_combined_R1.fq
gzip ./'PICO'$t'_R1.fq'
gunzip ./'PICO'$t'_R2.fq.gz'
cat ./'PICO'$t'_R2.fq' > ./PICO_combined_R2.fq
gzip ./'PICO'$t'_R2.fq'