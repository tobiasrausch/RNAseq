#!/bin/bash

export PATH=/g/funcgen/bin:${PATH}

# Download HG00096 FASTQ files
wget "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR188040/ERR188040_1.fastq.gz"
wget "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR188040/ERR188040_2.fastq.gz"

# Run RNA-seq processing
../src/align.sh HG00096 ERR188040_1.fastq.gz ERR188040_2.fastq.gz
