#!/bin/bash

if [ $# -ne 4 ]
then
    echo ""
    echo "Usage: $0 <read1.fq.gz> <read2.fq.gz> <output prefix> <strand [0|1|2]>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

# Activate environment
export PATH=${BASEDIR}/../mamba/bin:${PATH}

# CMD params
THREADS=4
FQ1=${1}
FQ2=${2}
OP=${3}
STRAND=${4}
HG=${BASEDIR}/../genome/GRCh38.fa

# Fastqc
mkdir -p ${OP}.read1.fastqc/
fastqc -t ${THREADS} -o ${OP}.read1.fastqc/ ${FQ1}
mkdir -p ${OP}.read2.fastqc/
fastqc -t ${THREADS} -o ${OP}.read2.fastqc/ ${FQ2}

# STAR alignment
STAR --runThreadN ${THREADS} --outFileNamePrefix ${OP}.star --outTmpDir ${OP}.tmpSTAR --genomeDir ${BASEDIR}/../genome/ --readFilesIn ${FQ1} ${FQ2} --readFilesCommand zcat
rm -rf ${OP}.tmpSTAR*

# Convert to BAM
samtools view -@ ${THREADS} -SbT ${HG} ${OP}.starAligned.out.sam > ${OP}.starraw.bam
rm ${OP}.starAligned.out.sam

# Sort & Index
samtools sort -@ ${THREADS} -o ${OP}.star.bam ${OP}.starraw.bam
samtools index ${OP}.star.bam
rm ${OP}.starraw.bam

# Subset to chr of interest
samtools view -b ${OP}.star.bam 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y > ${OP}.star.tmp.bam
mv ${OP}.star.tmp.bam ${OP}.star.bam
samtools index ${OP}.star.bam

# Basic alignment QC
samtools idxstats ${OP}.star.bam > ${OP}.idxstats
samtools flagstat ${OP}.star.bam > ${OP}.flagstat

# Run QC and gene counting using Alfred (check which strand option applies to your RNA-Seq kit)
alfred qc -r ${HG} -o ${OP}.alfred.tsv.gz ${OP}.star.bam
alfred count_rna -s ${STRAND} -g ${BASEDIR}/../gtf/Homo_sapiens.GRCh38.93.gtf.gz -o ${OP}.gene.count ${OP}.star.bam
alfred count_rna -s ${STRAND} -n fpkm -g ${BASEDIR}/../gtf/Homo_sapiens.GRCh38.93.gtf.gz -o ${OP}.gene.fpkm ${OP}.star.bam
alfred count_rna -s ${STRAND} -n fpkm_uq -g ${BASEDIR}/../gtf/Homo_sapiens.GRCh38.93.gtf.gz -o ${OP}.gene.fpkm_uq ${OP}.star.bam

# Create browser tracks
alfred tracks -o ${OP}.bedGraph.gz ${OP}.star.bam
igvtools totdf ${OP}.bedGraph.gz ${OP}.tdf hg19

