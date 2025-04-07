#!/bin/sh

#SBATCH --partition=ncpu
#SBATCH --job-name=cobalt
#SBATCH --mem=16G
#SBATCH -n 1
#SBATCH --time=168:00:00

# Script for processing cobaltseq samples
# find_files.sh should be run before this script

# Load modules
ml purge
ml Python/3.11.5-GCCcore-13.2.0

# Define inputs
PROJDIR=/nemo/stp/babs/working/bootj/projects/swantonc/eva.gongross/eg879
R1_FILES=${PROJDIR}/R1_files.txt
R2_FILES=${PROJDIR}/R2_files.txt
REPO=/nemo/stp/babs/working/bootj/github/cobaltSeq
BARCODES=${REPO}/barcodes.csv
#ITER=$(wc -l < ${R1_FILES})
ITER=5

# For loop to iterate over all samples
for i in $(seq ${ITER}); do

    # Define R1 and R2
    R1=$(sed -n "${i}p" ${R1_FILES})
    R2=$(sed -n "${i}p" ${R2_FILES})

    # Log
    echo "Iteration: ${i}, files:"
    echo ${R1}
    echo ${R2}

    # Run
    python ${REPO}/extract_reads_by_barcode.py ${R1} ${R2} ${BARCODES}

done