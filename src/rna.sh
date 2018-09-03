#!/bin/bash

if [ $# -ne 4 ]
then
    echo "**********************************************************************"
    echo "RNA-Seq analysis pipeline."
    echo "This program comes with ABSOLUTELY NO WARRANTY."
    echo ""
    echo "Version: 0.1.1"
    echo "Contact: Tobias Rausch (rausch@embl.de)"
    echo "**********************************************************************"
    echo ""
    echo "Usage: $0 <read1.fq.gz> <read2.fq.gz> <output prefix> <strand [0|1|2]>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

# CMD parameters
FQ1=${1}
FQ2=${2}
OUTP=${3}
STRAND=${4}
HG=${BASEDIR}/../genome/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa

# Align
${BASEDIR}/align.sh ${FQ1} ${FQ2} ${OUTP} ${STRAND}

# Call variants
${BASEDIR}/call.sh ${HG} ${OUTP}.star.bam ${OUTP}

# Phase variants
${BASEDIR}/phase.sh ${OUTP}.norm.filtered.vcf.gz ${OUTP}

# ASE using phased variants
${BASEDIR}/ase.sh ${HG} ${OUTP}.star.bam ${OUTP}.eagle2.bcf ${OUTP}
