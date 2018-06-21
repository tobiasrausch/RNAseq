# RNAseq
RNA-Seq processing pipeline


Prerequisites
-------------

Download STAR GRCh37.75 genome

`cd genome/ && ./starGenome.sh`

Download the GTF annotation file

`cd gtf/ && ./downloadGTF.sh`

Download the 1000 Genomes Reference Panel for allele-specific expression

`cd refpanel/ && ./download_1kGP_hg19.sh`


RNA-Seq Alignment and Gene Counting
-----------------------------------

`./src/rna.sh ERR188040_1.fastq.gz ERR188040_2.fastq.gz HG00096`

Example Workflow
----------------

To run an entire sample from the gEUVADIS data set

`cd example/ && ./gEUVADIS.sh`

