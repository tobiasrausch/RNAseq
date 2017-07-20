# RNAseq
RNA-Seq processing pipeline


Prerequisites
-------------

Download STAR GRCh37.75 genome

`cd genome/ && ./starGenome.sh`

Build a BED file of all exonic regions

`cd R && Rscript exon.R`

Download the GTF annotation file

`cd gtf/ && ./downloadGTF.sh`
