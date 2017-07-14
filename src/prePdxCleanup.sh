#!/bin/bash

if [ $# -ne 4 ]
then
    echo "**********************************************************************"
    echo "RNA-Seq processing pipeline"
    echo "This program comes with ABSOLUTELY NO WARRANTY."
    echo ""
    echo "Contact: Tobias Rausch (rausch@embl.de)"
    echo "**********************************************************************"
    echo ""
    echo "Usage: $0 <read1.fq.gz> <read2.fq.gz> <human_mouse.fa> <outprefix>"
    echo ""
    exit -1
fi


SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

export PATH=/g/funcgen/bin/:${PATH}

# CMD parameters
FQ1=${1}
FQ2=${2}
HG=${3}
OP=${4}
THREADS=4

# BWA
bwa mem -t ${THREADS} ${HG} ${FQ1} ${FQ2} | samtools view -@ ${THREADS} -bT ${HG} - > ${OP}.mouse.human.bam
samtools sort -@ ${THREADS} -o ${OP}.mouse.human.srt.bam ${OP}.mouse.human.bam
rm ${OP}.mouse.human.bam
samtools index ${OP}.mouse.human.srt.bam

# Mouse fraction
samtools idxstats ${OP}.mouse.human.srt.bam > ${OP}.mouse.human.idxstats
samtools flagstat ${OP}.mouse.human.srt.bam > ${OP}.mouse.human.flagstat

# Subset to human chr
samtools view -F 1804 -f 2 -b ${OP}.mouse.human.srt.bam `cat ${OP}.mouse.human.idxstats | grep "^chr" | cut -f 1 | tr '\n' ' '` > ${OP}.mouse.human.srt.clean.bam
samtools index ${OP}.mouse.human.srt.clean.bam
rm ${OP}.mouse.human.srt.bam ${OP}.mouse.human.srt.bam.bai

# Convert to FASTQ
samtools fastq -1 ${OP}.cleaned.1.fq -2 ${OP}.cleaned.2.fq ${OP}.mouse.human.srt.clean.bam
rm ${OP}.mouse.human.srt.clean.bam ${OP}.mouse.human.srt.clean.bam.bai
gzip ${OP}.cleaned.1.fq ${OP}.cleaned.2.fq

# Run main RNA-seq pipeline
${BASEDIR}/rna.sh ${OP}.cleaned.1.fq.gz ${OP}.cleaned.2.fq.gz ${OP}
