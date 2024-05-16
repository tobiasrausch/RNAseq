#!/bin/bash

if [ ! -f Homo_sapiens.GRCh38.93.gtf.gz ]
then
    wget 'ftp://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz'
fi
zcat Homo_sapiens.GRCh38.93.gtf.gz | awk '($2=="havana" || $2=="ensembl_havana") && $3=="exon"' | cut -f 1,4,5 | grep "^[0-9X]" | sort -k1,1V -k2,2n | uniq > exons.bed
bedtools merge -i exons.bed > coding.hg38.bed 
cat coding.hg38.bed | awk '{SUM+=($3-$2);} END {print SUM;}'
rm exons.bed

