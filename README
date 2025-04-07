# README

## Sequencing folder
`/nemo/stp/babs/inputs/sequencing/data/swantonc/su-kit.chew/DN20060/primary_data/`

## Sequencing runs
- `201008_K00371_0408_AHHYFWBBXY`
- `210831_A01366_0049_AH35HJDMXY`
- `210902_A01366_0050_AH2WG2DMXY`
- `211116_A01366_0091_BH3YLGDMXY`
- `211130_A01366_0101_BH3YWKDMXY`
- `220125_A01366_0132_AH5VMHDMXY`
- `220422_A01366_0181_AHCH3FDMXY`
- `220428_A01366_0186_AHCCVFDMXY`
- `220509_A01366_0195_AHCCTWDMXY`
- `220610_A01366_0217_BHCGJGDMXY`
- `230428_A01366_0383_BH3HYNDSX7`
- `230505_A01366_0385_BHYC5YDSX5`

## extract_reads_by_barcode.py

### Overview:
The script will go through the paired `fastq` files (R1 & R2) and extract reads containing the all the gRNA barcodes,  supplied in `<barcodes.csv>`, and save reads to individual files for each gRNA. It will also extract and count occurences of the random 30nt barcode that immediately follows the gRNA barcode. the random barcode must start and end with `GC`, and also contain `TA` at positions 7 and 21, and another `GC` at position 14. Output files will be pre-appended with the barcode/sgRNA name.

NOTE: input files must be the full absolute path to the file, the script makes assumptions on NEMO directory structure.

### Usage: 
`python extract_reads_by_barcode.py <R1.fastq(.gz)> <R2.fastq(.gz)> <barcodes.csv>`

### Details:
`R1.fastq(.gz)` is the full absolute path to a R1 fastq containing raw reads

`R2.fastq(.gz)` is the full absolute path to a R2 fastq containing raw reads

`barcodes.csv` should contain two columns, col1 should be the barcode/sgRNA name, col2 should be the actual barcode

## find_files.sh

Bash script for finding and sorting all R1 and R2 files ready for processing by sbatch script. The search directory (`SEARCHDIR`) should be specified in the script. Output will be 2 sorted files listing all R1 and R2 files.

## sbatch.sh

Bash script for submitting all files for processing by the `extract_reads_by_barcode.py` script. `find_files.sh` should be run before running this script.

If running testing set the `ITER` variable to the number of files you want to test on. The top `n` files will be used from the top of the R1 and R2 list files.

## Is R2 just reverse complement of R1?

It doesn't look like it...

Search for reverse complement of barcodes in fastq file

9p21_01 = CCACTTCT

9p21_01 rev complement = AGAAGTGG
