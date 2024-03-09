#!/bin/bash

# Переходим в текущую директорию
cd "$(dirname "$0")"

# Проходимся по всем директориям в текущей
for DIR in */; do
  if [ -d "$DIR" ]; then
    # Удаляем все директории, кроме fastq
    find "$DIR" -mindepth 1 -maxdepth 1 -type d -not -name "fastq" -exec rm -rf {} \;

    # Запускаем команду в каждой директории
    bash "/home/aapopov/tool/juicer1.6_compact-main/scripts/juicer.sh" -g hg38 -d /data/aapopov/Projects/HumanDros/Zapusk_Sirius_Human_Dros/${DIR} -s none -p "/data/aapopov/Projects/HumanDros/CTCF_Aux_data_for_juicer/genome/Homo_sapiens.GRCh38.dna.primary_assembly.chrom.sizes" -z "/data/aapopov/Projects/HumanDros/CTCF_Aux_data_for_juicer/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa" -D "/home/aapopov/tool/juicer1.6_compact-main" -t 20
  fi
done