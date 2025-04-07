#!/bin/sh

#SBATCH --partition=ncpu
#SBATCH --job-name=findfiles
#SBATCH --mem=4G
#SBATCH -n 1
#SBATCH --time=1:00:00
#SBATCH --error=findfiles.err

SEARCHDIR=/nemo/stp/babs/inputs/sequencing/data/swantonc/su-kit.chew/DN20060/primary_data

find -L ${SEARCHDIR} -name "*R1_001.fastq.gz" | sort >> R1_files.txt
find -L ${SEARCHDIR} -name "*R2_001.fastq.gz" | sort >> R2_files.txt