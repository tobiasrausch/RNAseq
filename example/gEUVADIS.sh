#!/bin/bash

# Download HG00096 FASTQ files
wget "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR188040/ERR188040_1.fastq.gz"
wget "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR188040/ERR188040_2.fastq.gz"

# Run RNA-seq processing
../src/rna.sh ERR188040_1.fastq.gz ERR188040_2.fastq.gz HG00096
