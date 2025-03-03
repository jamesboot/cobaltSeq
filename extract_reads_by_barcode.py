import gzip
import csv
import os
import sys

def read_fastq(file_path):
    """
    Generator to read a FASTQ file (can be gzipped).
    Yields tuples: (header, sequence, plus, quality).
    """
    if file_path.endswith('.gz'):
        opener = gzip.open
        mode = 'rt'
    else:
        opener = open
        mode = 'r'

    with opener(file_path, mode) as f:
        while True:
            header = f.readline().strip()
            sequence = f.readline().strip()
            plus = f.readline().strip()
            quality = f.readline().strip()
            if not header:
                break
            yield header, sequence, plus, quality

def load_barcodes(csv_file):
    """
    Load barcodes from a CSV file.
    Returns a dictionary {barcode_name: barcode_sequence}
    """
    barcodes = {}
    with open(csv_file, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) != 2:
                continue
            name, seq = row
            barcodes[name] = seq
    return barcodes

def extract_reads_for_barcodes(input_fastq, barcodes_csv):
    """
    Extract reads containing each barcode from the CSV file.
    Each barcode gets its own FASTQ file.
    """
    # Load barcodes from the CSV
    barcodes = load_barcodes(barcodes_csv)

    # Determine the base name of the input file (remove .gz/.fastq extensions)
    base_name = os.path.basename(input_fastq)
    if base_name.endswith('.gz'):
        base_name = base_name[:-3]
    if base_name.endswith('.fastq'):
        base_name = base_name[:-6]

    # Prepare output file handles for each barcode
    output_handles = {}
    for barcode_name in barcodes:
        output_file = f'{base_name}_{barcode_name}.fastq'
        output_handles[barcode_name] = open(output_file, 'w')

    # Scan through FASTQ file and write reads to matching barcode files
    for header, sequence, plus, quality in read_fastq(input_fastq):
        for barcode_name, barcode_seq in barcodes.items():
            if barcode_seq in sequence:
                out_f = output_handles[barcode_name]
                out_f.write(f'{header}\n{sequence}\n{plus}\n{quality}\n')
                break  # Stop after the first matching barcode (optional, depends if you want multi-matches)

    # Close all output files
    for handle in output_handles.values():
        handle.close()

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python extract_reads_by_barcode.py <input.fastq(.gz)> <barcodes.csv>")
        sys.exit(1)

    input_fastq = sys.argv[1]
    barcodes_csv = sys.argv[2]

    extract_reads_for_barcodes(input_fastq, barcodes_csv)