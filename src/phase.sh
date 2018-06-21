#!/bin/bash

if [ $# -ne 2 ]
then
    echo ""
    echo "Usage: $0 <input.vcf.gz> <outprefix>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

# Activate environment
export PATH=${BASEDIR}/../bin/Eagle_v2.4:${BASEDIR}/../bin/bin:${PATH}
source activate ${BASEDIR}/../bin/envs/rna

# CMD params
THREADS=4
VCF=${1}
OP=${2}

# Fetch input variants
bcftools view -m2 -M2 --types snps,indels ${VCF} | sed 's/^##contig=<ID=chr/##contig=<ID=/' | sed 's/^chr//' | bcftools annotate -O b -o ${OP}.input.bcf -x INFO,^FORMAT/GT -
bcftools index ${OP}.input.bcf

# Phase against panel
FILES=""
#for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
for CHR in 1 2
do
    echo "Eagle2 phasing chr${CHR}"
    if [ `bcftools view ${OP}.input.bcf ${CHR} | grep -m 1 "^#CHROM" -A 1 | wc -l` -eq 2 ]
    then
	eagle --numThreads ${THREADS} --vcfRef ${BASEDIR}/../refpanel/chr${CHR}.bcf --vcfTarget ${OP}.input.bcf --geneticMapFile ${BASEDIR}/../refpanel/genetic_map_hg19_withX.txt.gz --outPrefix ${OP}.chr${CHR}.eagle2 --vcfOutFormat b --chrom ${CHR} 2>&1 | gzip -c > ${OP}.chr${CHR}.eagle2.log.gz
	bcftools index ${OP}.chr${CHR}.eagle2.bcf
	FILES=${FILES}" "${OP}.chr${CHR}.eagle2.bcf
    fi
done
rm ${OP}.input.bcf ${OP}.input.bcf.csi
bcftools concat -O b -o ${OP}.eagle2.bcf ${FILES}
bcftools index ${OP}.eagle2.bcf

# Clean-up
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
    if [ -f ${OP}.chr${CHR}.eagle2.bcf ]
    then
	rm ${OP}.chr${CHR}.eagle2.bcf ${OP}.chr${CHR}.eagle2.bcf.csi ${OP}.chr${CHR}.eagle2.log.gz
    fi
done

# Deactivate environment
source deactivate
