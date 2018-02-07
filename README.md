# RNAseq
RNA-Seq processing pipeline


Prerequisites
-------------

Download STAR GRCh37.75 genome

`cd genome/ && ./starGenome.sh`

Download the GTF annotation file

`cd gtf/ && ./downloadGTF.sh`



RNA-Seq Alignment and Gene Counting
-----------------------------------

`./src/align.sh HG00096 ERR188040_1.fastq.gz ERR188040_2.fastq.gz`


Example Workflow
----------------

Align one sample from the gEUVADIS data set

`cd example/ && ./gEUVADIS.sh`

