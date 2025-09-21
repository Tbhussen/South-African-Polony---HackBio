#!/bin/bash

#set -e

first_part="SRR270133"
middle_part="_Genome_Sequencing_of_Listeria_monocytogenes_SA_outbreak_2017_"
last_part_1="1.fastq.gz"
last_part_2="2.fastq.gz"

mkdir -p qc
mkdir -p trimmed
mkdir -p assembly
mkdir -p AMR

# Run FastQC on all FASTQ files in current dir
for file in raw/*.fastq.gz; do
    fastqc "$file" -o qc/
done

# Run fastp on samples SRR27013311 ... SRR27013343
for i in {17..43}; do
    fastp \
    -i raw/${first_part}${i}${middle_part}${last_part_1} \
    -I raw/${first_part}${i}${middle_part}${last_part_2} \
    -o trimmed/${first_part}${i}${middle_part}.trim.${last_part_1} \
    -O trimmed/${first_part}${i}${middle_part}.trim.${last_part_2}

    if [ $? -ne 0 ]; then
        echo "Fastp failed for $sample, skipping..."
        rm -rf assembly/"$sample"
        continue
    fi

done

# Assemble each pair + Detect AMR genes
for fq in trimmed/*_.trim.1.fastq.gz; do
    s=$(basename "$fq" _.trim.1.fastq.gz)
    spades.py \
      -1 trimmed/${s}_.trim.1.fastq.gz \
      -2 trimmed/${s}_.trim.2.fastq.gz \
      -o assembly/${s} || true

    quast.py \
      assembly/${s}/contigs.fasta -o quast_report

    abricate \
      assembly/${s}/contigs.fasta \
      > AMR/resistance_report.txt
done

# Apply abricate with VFDB to detect toxins on a single read
abricate --db vfdb assembly/SRR27013312_Genome_Sequencing_of_Listeria_monocytogenes_SA_outbreak_2017/contigs.fasta > VF/SRR27013312_Genome_Sequencing_of_Listeria_monocytogenes_SA_outbreak_2017_vfdb.tab
