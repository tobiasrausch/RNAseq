#!/bin/bash

export PATH=/g/funcgen/bin:${PATH}

if [ $# -ne 1 ]
then
    echo "**********************************************************************"
    echo "RNA-Seq processing pipeline."
    echo "This program comes with ABSOLUTELY NO WARRANTY."
    echo ""
    echo "Contact: Tobias Rausch (rausch@embl.de)"
    echo "**********************************************************************"
    echo ""
    echo "Usage: $0 <hg19.fa>"
    echo ""
    exit -1
fi
    

HG=${1}

# Download HG00096 FASTQ files
wget "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR188040/ERR188040_1.fastq.gz"
wget "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR188040/ERR188040_2.fastq.gz"

# Run RNA-seq processing
../src/rna.sh ERR188040_1.fastq.gz ERR188040_2.fastq.gz HG00096
