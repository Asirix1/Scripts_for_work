#!/bin/bash

# Заходим во все поддиректории
for dir in $(find . -type d); do
    # Проверяем наличие файлов *R1.fastq и *R2.fastq
    if [ -f "$dir"/*R1.fastq ] && [ -f "$dir"/*R2.fastq ]; then
        # Проходим через каждые четыре строки и добавляем \1 и \2 к первому столбцу
        for file in "$dir"/*R1.fastq; do
            awk '{if(NR % 4 == 1) {print $1"\\1 "$2} else {print $0}}' "$file" > temp.fastq
            mv temp.fastq "$file"
        done
        for file in "$dir"/*R2.fastq; do
            awk '{if(NR % 4 == 1) {print $1"\\2 "$2} else {print $0}}' "$file" > temp.fastq
            mv temp.fastq "$file"
        done
    fi
done