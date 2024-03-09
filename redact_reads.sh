#!/bin/bash

# ������� �� ��� �������������
for dir in $(find . -type d); do
    # ��������� ������� ������ *R1.fastq � *R2.fastq
    if [ -f "$dir"/*R1.fastq ] && [ -f "$dir"/*R2.fastq ]; then
        # �������� ����� ������ ������ ������ � ��������� \1 � \2 � ������� �������
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