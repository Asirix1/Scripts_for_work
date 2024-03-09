#!/bin/bash

# ��������� � ������� ����������
cd "$(dirname "$0")"

# ���� ��� ����������
for dir in */; do
    if [ -d "$dir" ]; then
        cd "$dir"
        
        # ��������� ������� ���������� fastq
        if [ -d "fastq" ]; then
            cd "fastq"
            
            # ���������� ����� R1.gz � ����
            cat *R1* > "${dir%/}_R1.fastq"
            gunzip "${dir%/}_R1.fastq.gz"
            # ���������� ����� R2.gz � ����
            cat *R2* > "${dir%/}_R2.fastq"
            gunzip "${dir%/}_R2.fastq.gz"
            cd ..
        fi
        cd ..
    fi
done
