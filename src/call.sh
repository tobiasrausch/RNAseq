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
export PATH=${BASEDIR}/../bin/bin:${PATH}
source activate ${BASEDIR}/../bin/envs/rna

# CMD params
THREADS=4
GENOME=${1}
BAM=${2}
OUTP=${3}

# Freebayes
freebayes --no-indels --no-mnps --no-complex --no-partial-observations --min-repeat-entropy 1 --report-genotype-likelihood-max --min-alternate-fraction 0.15 --fasta-reference ${GENOME} --genotype-qualities ${BAM} -v ${OUTP}.vcf
bgzip ${OUTP}.vcf
tabix ${OUTP}.vcf.gz

# Normalize VCF
bcftools norm -O z -o ${OUTP}.norm.vcf.gz -f ${GENOME} -m -both ${OUTP}.vcf.gz
tabix ${OUTP}.norm.vcf.gz
rm ${OUTP}.vcf.gz ${OUTP}.vcf.gz.tbi

# Fixed threshold filtering
bcftools filter -O z -o ${OUTP}.norm.filtered.vcf.gz -e '%QUAL<=20 || %QUAL/AO<=2 || SAF<=2 || SAR<=2' ${OUTP}.norm.vcf.gz
tabix ${OUTP}.norm.filtered.vcf.gz
rm ${OUTP}.norm.vcf.gz ${OUTP}.norm.vcf.gz.tbi

# Deactivate environment
source deactivate
