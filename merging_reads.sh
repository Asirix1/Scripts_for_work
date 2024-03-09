#!/bin/bash

# Переходим в текущую директорию
cd "$(dirname "$0")"

# Ищем все директории
for dir in */; do
    if [ -d "$dir" ]; then
        cd "$dir"
        
        # Проверяем наличие директории fastq
        if [ -d "fastq" ]; then
            cd "fastq"
            
            # Объединяем файлы R1.gz в один
            cat *R1* > "${dir%/}_R1.fastq"
            gunzip "${dir%/}_R1.fastq.gz"
            # Объединяем файлы R2.gz в один
            cat *R2* > "${dir%/}_R2.fastq"
            gunzip "${dir%/}_R2.fastq.gz"
            cd ..
        fi
        cd ..
    fi
done
