import gzip
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path

# Function for count total reads in a FASTQ file
# This function counts the number of reads in a FASTQ file by counting the number of lines
# and dividing by 4 (since each read consists of 4 lines: header, sequence, plus line, quality).
# It handles both gzipped and uncompressed files.
def count_total_reads(fastq_file):
    """Count the total number of reads in a FASTQ file."""
    if fastq_file.endswith('.gz'):
        opener = gzip.open
        mode = 'rt'
    else:
        opener = open
        mode = 'r'

    count = 0
    with opener(fastq_file, mode) as f:
        for _ in f:
            count += 1
    return count // 4  # Each read consists of 4 lines

# Function to read paired-end FASTQ files
# This function reads paired-end FASTQ files (R1 and R2) line by line.
# It yields a tuple containing the header, sequence, plus line, and quality line for both R1 and R2.
# It handles both gzipped and uncompressed files.
# The function uses a generator to yield the reads one by one, which is memory efficient for large files.
# It assumes that the input files are in the correct format and that R1 and R2 are in the same order.
# The function will stop reading when it reaches the end of either file.
def read_paired_fastq(r1_path, r2_path):
    """Generator to read paired-end FASTQ files (R1 and R2)."""
    # Check if the input files are gzipped or not
    # and set the appropriate opener and mode
    if r1_path.endswith('.gz'):
        opener = gzip.open
        mode = 'rt'
    else:
        opener = open
        mode = 'r'
        
    # Open both files using the appropriate opener
    # and read them line by line
    with opener(r1_path, mode) as f1, opener(r2_path, mode) as f2:
        while True:
            r1_header = f1.readline().strip()
            r1_sequence = f1.readline().strip()
            r1_plus = f1.readline().strip()
            r1_quality = f1.readline().strip()

            r2_header = f2.readline().strip()
            r2_sequence = f2.readline().strip()
            r2_plus = f2.readline().strip()
            r2_quality = f2.readline().strip()

            if not r1_header or not r2_header:
                break

            yield (r1_header, r1_sequence, r1_plus, r1_quality, 
                   r2_header, r2_sequence, r2_plus, r2_quality)

# Function to load barcodes from a CSV file
def load_barcodes(csv_file):
    """Load barcodes from a CSV file."""
    barcodes = {}
    with open(csv_file, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) != 2:
                continue
            name, seq = row
            barcodes[name] = seq
    return barcodes

# Main function to extract paired reads and 30-nt sequences starting with 'GC'
# This function takes the input R1 and R2 FASTQ files, a CSV file with barcodes,
# and extracts the paired reads and 30-nt sequences starting with 'GC'.
# It creates output files for each barcode and writes the corresponding reads and sequences.
# It also counts the total number of reads in the input files and writes the counts to a CSV file.
# The output files are named based on the sample ID, Illumina ID, lane ID, and barcode name.
def extract_paired_reads_and_gc_30nt(input_r1, input_r2, barcodes_csv):
    """Extract paired reads and 30-nt sequences starting with 'GC'."""
    barcodes = load_barcodes(barcodes_csv)

    # Define elements of sample name 
    name_elements = os.path.basename(input_r1).split('_')
    sample_id = name_elements[0]
    illumina_id = name_elements[1]
    lane_id = name_elements[2]
    run_id = input_r1.split('/')[-3] # Assumes NEMO structure - run folder is here

    # Create output folder 
    Path(run_id).mkdir(parents=True, exist_ok=True)

    # Initialize dictionaries to store output handles and sequence counts
    output_handles = {}
    sequence_files = {}
    sequence_counts = defaultdict(lambda: defaultdict(int))
    gRNA_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))

    # Create output files for each barcode
    for barcode_name in barcodes:
        r1_output_file = f'{run_id}/{sample_id}_{illumina_id}_{lane_id}_{barcode_name}_R1.fastq'
        r2_output_file = f'{run_id}/{sample_id}_{illumina_id}_{lane_id}_{barcode_name}_R2.fastq'
        seq_output_file = f'{run_id}/{sample_id}_{illumina_id}_{lane_id}_{barcode_name}_30nt.txt'

        output_handles[barcode_name] = (open(r1_output_file, 'w'), open(r2_output_file, 'w'))
        sequence_files[barcode_name] = open(seq_output_file, 'w')

    # Read paired-end FASTQ files and extract reads and sequences
    # for each barcode
    for (r1_header, r1_sequence, r1_plus, r1_quality, 
         r2_header, r2_sequence, r2_plus, r2_quality) in read_paired_fastq(input_r1, input_r2):

        # Check if the read contains any of the barcodes
        # and extract the corresponding reads and sequences
        for barcode_name, barcode_seq in barcodes.items():
            # Check if the barcode sequence is present in the R1 sequence
            # and find its index
            index = r1_sequence.find(barcode_seq)
            # If the barcode sequence is found, extract the reads
            # and sequences and write them to the corresponding files
            if index != -1:
                r1_out, r2_out = output_handles[barcode_name]
                seq_out = sequence_files[barcode_name]

                r1_out.write(f'{r1_header}\n{r1_sequence}\n{r1_plus}\n{r1_quality}\n')
                r2_out.write(f'{r2_header}\n{r2_sequence}\n{r2_plus}\n{r2_quality}\n')
                
                # Increment the count of valid barcodes
                gRNA_counts[sample_id][barcode_name]['valid_gRNA'] += 1

                # Extract the 30-nt sequence starting from the position
                # after the barcode sequence
                # Check if the next 30-nt sequence starts with 'GC'
                # and contains the specified patterns
                # and ends with 'GC'
                start_pos = index + len(barcode_seq)
                if start_pos + 30 <= len(r1_sequence):
                    next_30nt = r1_sequence[start_pos:start_pos + 30]
                    if next_30nt.startswith("GC") and next_30nt[7:9] == 'TA' and next_30nt[14:16] == 'GC' and next_30nt[21:23] == 'TA' and next_30nt.endswith("GC"):
                        # Write the 30-nt sequence to the corresponding file
                        # and update the count
                        seq_out.write(next_30nt + "\n")
                        sequence_counts[barcode_name][next_30nt] += 1
                        gRNA_counts[sample_id][barcode_name]['valid_rand_bc'] += 1
                    else:
                        # If the sequence does not match the criteria, increment the invalid count
                        gRNA_counts[sample_id][barcode_name]['rand_bc_seq_error'] += 1
                else:
                    # If the sequence is too short, increment the invalid count
                    gRNA_counts[sample_id][barcode_name]['rand_bc_too_short'] += 1    
                # Break because we only want to process one barcode per read
                # If the barcode sequence is found, we can stop checking other barcodes
                break
            
    # Close all output files            
    for r1_out, r2_out in output_handles.values():
        r1_out.close()
        r2_out.close()
    for seq_out in sequence_files.values():
        seq_out.close()
        
    # Count total reads in input FASTQ files
    total_reads_r1 = count_total_reads(input_r1)
    total_reads_r2 = count_total_reads(input_r2)

    # Write the counts to a CSV file for each barcode
    # The CSV file contains the sequence, count, barcode name, sample ID,
    # Illumina ID, lane ID, run ID, and total reads in R1 and R2
    # The output files are named based on the sample ID, Illumina ID,
    # lane ID, and barcode name
    for barcode_name in sequence_counts:
        count_file = f"{run_id}/{sample_id}_{illumina_id}_{lane_id}_{barcode_name}_counts.csv"
        with open(count_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            # Write the header for the CSV file
            writer.writerow(['rand_bc', 'count', 'gRNA', 'sample_id', 'illumina_id', 'lane_id', 'run_id', 'total_reads_R1', 'total_reads_R2'])
            # Write the counts for each sequence in descending order
            # sorted by count
            for sequence, count in sorted(sequence_counts[barcode_name].items(), key=lambda x: -x[1]):
                writer.writerow([sequence, count, barcode_name, sample_id, illumina_id, lane_id, run_id, total_reads_r1, total_reads_r2])  # No header, includes barcode name

    # Write the gRNA counts to a CSV file
    gRNA_count_file = f"{run_id}/{sample_id}_{illumina_id}_{lane_id}_gRNA_counts.csv"
    with open(gRNA_count_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Write the header for the gRNA counts CSV file
        writer.writerow(['sample_id', 'barcode_name', 'total_reads_R1', 'total_reads_R2', 'valid_gRNA', 'valid_rand_bc', 'rand_bc_seq_error', 'rand_bc_too_short'])
        # Write the counts for each barcode
        for barcode_name, counts in gRNA_counts[sample_id].items():
            writer.writerow([sample_id, barcode_name, total_reads_r1, total_reads_r2, counts['valid_gRNA'], counts['valid_rand_bc'], counts['rand_bc_seq_error'], counts['rand_bc_too_short']])

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python extract_reads_by_barcode.py <R1.fastq(.gz)> <R2.fastq(.gz)> <barcodes.csv>")
        sys.exit(1)

    input_r1 = sys.argv[1]
    input_r2 = sys.argv[2]
    barcodes_csv = sys.argv[3]

    extract_paired_reads_and_gc_30nt(input_r1, input_r2, barcodes_csv)
