#!/bin/sh

#SBATCH --partition=ncpu
#SBATCH --job-name=concatfiles
#SBATCH --mem=4G
#SBATCH -n 1
#SBATCH --time=1:00:00

# Script to concatenate output count files outut from extract_reads_by_barcode.py

# Define directory to search
SEARCHDIR=/nemo/stp/babs/working/bootj/projects/swantonc/eva.gongross/eg879

# Define output file
OUTPUT_FILE1="allcounts.csv"
OUTPUT_FILE2="all_gRNA_counts.csv"

# Find files 
find -L ${SEARCHDIR} -name "*_counts.csv" | sort >> all_counts_outs.txt
find -L ${SEARCHDIR} -name "*_gRNA_counts.csv" | sort >> all_gRNA_counts_outs.txt

# Concatenate the counts files
# Create or clear the output file
> "${OUTPUT_FILE1}"
# Read each file from the list and append it to the output file
# Skip first line of each file
while IFS= read -r file; do
    if [ -f "$file" ]; then
        # Skip the header line of each file except the first one
        if [ "$file" != "$(head -n 1 all_counts_outs.txt)" ]; then
            tail -n +2 "$file" >> "${OUTPUT_FILE1}"
        else
            cat "$file" >> "${OUTPUT_FILE1}"
        fi
    else
        echo "Warning: File '$file' not found, skipping."
    fi
done < all_counts_outs.txt

# Concatenate the gRNA counts files
# Create or clear the output file
> "${OUTPUT_FILE2}"
# Read each file from the list and append it to the output file
# Skip first line of each file
while IFS= read -r file; do
    if [ -f "$file" ]; then
        # Skip the header line of each file except the first one
        if [ "$file" != "$(head -n 1 all_gRNA_counts_outs.txt)" ]; then
            tail -n +2 "$file" >> "${OUTPUT_FILE2}"
        else
            cat "$file" >> "${OUTPUT_FILE2}"
        fi
    else
        echo "Warning: File '$file' not found, skipping."
    fi
done < all_gRNA_counts_outs.txt
