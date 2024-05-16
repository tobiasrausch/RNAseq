#!/bin/bash

if [ $# -ne 3 ]
then
    echo ""
    echo "Usage: $0 <genome.fa> <aligned.bam> <output prefix>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

# Activate environment
export PATH=${BASEDIR}/../mamba/bin:${PATH}

# CMD params
THREADS=4
GENOME=${1}
BAM=${2}
OUTP=${3}

# Strelka2
python2.7 /opt/dev/RNASeq/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py --rna --referenceFasta=${GENOME} --bam=${BAM} --runDir=strelka_${OUTP}
python2.7 strelka_${OUTP}/runWorkflow.py -m local -j ${THREADS}

# Normalize VCF
bcftools view -f 'PASS,.' strelka_${OUTP}/results/variants/genome.S1.vcf.gz | bcftools norm -O b -o ${OUTP}.norm.bcf -f ${GENOME} -m -both -
bcftools index ${OUTP}.norm.bcf

# Subset to coding regions
bcftools view -T ${BASEDIR}/../gtf/coding.hg38.bed ${OUTP}.norm.bcf | bgzip > ${OUTP}.vcf.gz
tabix ${OUTP}.vcf.gz
rm ${OUTP}.norm.bcf ${OUTP}.norm.bcf.csi
rm -rf strelka_${OUTP}/
