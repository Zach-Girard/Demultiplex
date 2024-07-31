#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=compute
#SBATCH --cpus-per-task=10
#SBATCH --mem=100GB
#SBATCH --job-name=CountN_Index


/usr/bin/time -v zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | awk 'NR == 2 || NR % 4 == 2' | grep "N" | wc -l 
/usr/bin/time -v zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | awk 'NR == 2 || NR % 4 == 2' | grep "N" | wc -l 