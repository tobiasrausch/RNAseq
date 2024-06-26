#!/bin/bash

if [ $# -ne 3 ]
then
    echo "**********************************************************************"
    echo "RNA-Seq analysis pipeline."
    echo "This program comes with ABSOLUTELY NO WARRANTY."
    echo ""
    echo "Version: 0.1.1"
    echo "Contact: Tobias Rausch (rausch@embl.de)"
    echo "**********************************************************************"
    echo ""
    echo "Usage: $0 <read.fq.gz> <output prefix> <strand [0|1|2]>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

# CMD parameters
FQ=${1}
OUTP=${2}
STRAND=${3}
HG=${BASEDIR}/../genome/GRCh38.fa

# Align
if [ ! -f ${OUTP}.bam ]
then
    ${BASEDIR}/alignSE.sh ${FQ} ${OUTP} ${STRAND}
fi

# Call variants
if [ ! -f ${OUTP}.vcf.gz ]
then
    ${BASEDIR}/call.sh ${HG} ${OUTP}.star.bam ${OUTP}
fi

# Phase variants
${BASEDIR}/phase.sh ${OUTP}.norm.filtered.vcf.gz ${OUTP}

# ASE using phased variants
${BASEDIR}/ase.sh ${HG} ${OUTP}.star.bam ${OUTP}.eagle2.bcf ${OUTP}
