#!/bin/bash

sra_ids="SRR4068197"

for sra_id in ${sra_ids}; do
    fastq-dump --gzip --outdir fastq --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip ${sra_id}
done