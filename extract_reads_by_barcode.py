import gzip
import csv
import os
import sys
from collections import defaultdict
from pathlib import Path

def read_paired_fastq(r1_path, r2_path):
    """Generator to read paired-end FASTQ files (R1 and R2)."""
    if r1_path.endswith('.gz'):
        opener = gzip.open
        mode = 'rt'
    else:
        opener = open
        mode = 'r'

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

    output_handles = {}
    sequence_files = {}
    sequence_counts = defaultdict(lambda: defaultdict(int))

    for barcode_name in barcodes:
        r1_output_file = f'{run_id}/{sample_id}_{illumina_id}_{lane_id}_{barcode_name}_R1.fastq'
        r2_output_file = f'{run_id}/{sample_id}_{illumina_id}_{lane_id}_{barcode_name}_R2.fastq'
        seq_output_file = f'{run_id}/{sample_id}_{illumina_id}_{lane_id}_{barcode_name}_30nt.txt'

        output_handles[barcode_name] = (open(r1_output_file, 'w'), open(r2_output_file, 'w'))
        sequence_files[barcode_name] = open(seq_output_file, 'w')

    for (r1_header, r1_sequence, r1_plus, r1_quality, 
         r2_header, r2_sequence, r2_plus, r2_quality) in read_paired_fastq(input_r1, input_r2):

        for barcode_name, barcode_seq in barcodes.items():
            index = r1_sequence.find(barcode_seq)
            if index != -1:
                r1_out, r2_out = output_handles[barcode_name]
                seq_out = sequence_files[barcode_name]

                r1_out.write(f'{r1_header}\n{r1_sequence}\n{r1_plus}\n{r1_quality}\n')
                r2_out.write(f'{r2_header}\n{r2_sequence}\n{r2_plus}\n{r2_quality}\n')

                start_pos = index + len(barcode_seq)
                if start_pos + 30 <= len(r1_sequence):
                    next_30nt = r1_sequence[start_pos:start_pos + 30]
                    if next_30nt.startswith("GC") and next_30nt[7:9] == 'TA' and next_30nt[14:16] == 'GC' and next_30nt[21:23] == 'TA' and next_30nt.endswith("GC"):
                        # Write the 30-nt sequence to the corresponding file
                        # and update the count
                        seq_out.write(next_30nt + "\n")
                        sequence_counts[barcode_name][next_30nt] += 1
                
                break

    for r1_out, r2_out in output_handles.values():
        r1_out.close()
        r2_out.close()
    for seq_out in sequence_files.values():
        seq_out.close()

    for barcode_name in sequence_counts:
        count_file = f"{run_id}/{sample_id}_{illumina_id}_{lane_id}_{barcode_name}_counts.csv"
        with open(count_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            for sequence, count in sorted(sequence_counts[barcode_name].items(), key=lambda x: -x[1]):
                writer.writerow([sequence, count, barcode_name, sample_id, illumina_id, lane_id, run_id])  # No header, includes barcode name

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python extract_reads_by_barcode.py <R1.fastq(.gz)> <R2.fastq(.gz)> <barcodes.csv>")
        sys.exit(1)

    input_r1 = sys.argv[1]
    input_r2 = sys.argv[2]
    barcodes_csv = sys.argv[3]

    extract_paired_reads_and_gc_30nt(input_r1, input_r2, barcodes_csv)
