#!/bin/bash

if [ $# -lt 2 ]
then
    echo "**********************************************************************"
    echo "RNA-Seq processing pipeline"
    echo "This program comes with ABSOLUTELY NO WARRANTY."
    echo ""
    echo "Contact: Tobias Rausch (rausch@embl.de)"
    echo "**********************************************************************"
    echo ""
    echo "Usage: $0 <outprefix> <read1.fq.gz> [<read2.fq.gz>]"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

export PATH=/g/funcgen/bin/:${PATH}

# CMD params
THREADS=4
OP=${1}
FQ1=${2}
PE=0
if [ $# -eq 3 ]
then
    PE=1
    FQ2=${3}
fi
HG=${BASEDIR}/../genome/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa

# Fastqc
mkdir -p ${OP}.read1.fastqc/
fastqc -t ${THREADS} -o ${OP}.read1.fastqc/ ${FQ1}
if [ ${PE} -eq 1 ]
then
    mkdir -p ${OP}.read2.fastqc/
    fastqc -t ${THREADS} -o ${OP}.read2.fastqc/ ${FQ2}
fi

# STAR alignment
if [ ${PE} -eq 1 ]
then
    STAR --runThreadN ${THREADS} --outFileNamePrefix ${OP}.star --outTmpDir ${OP}.tmpSTAR --genomeDir ${BASEDIR}/../genome/ --readFilesIn ${FQ1} ${FQ2} --readFilesCommand zcat
else
    STAR --runThreadN ${THREADS} --outFileNamePrefix ${OP}.star --outTmpDir ${OP}.tmpSTAR --genomeDir ${BASEDIR}/../genome/ --readFilesIn ${FQ1} --readFilesCommand zcat
fi
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

# Run QC and gene counting using Alfred
alfred qc -r ${HG} -o ${OP}.alfred.tsv.gz ${OP}.star.bam
alfred count_rna -g ${BASEDIR}/../gtf/Homo_sapiens.GRCh37.75.gtf.gz -o ${OP}.gene.count ${OP}.star.bam

# Fix chromosome names
samtools view -H ${OP}.star.bam > ${OP}.header
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
do
    sed -i "s/SN:${CHR}\t/SN:chr${CHR}\t/" ${OP}.header
done
cat ${OP}.header | grep -v -P "^@SQ\tSN:[A-Z]" > ${OP}.header.tmp
mv ${OP}.header.tmp ${OP}.header
samtools reheader ${OP}.header ${OP}.star.bam > ${OP}.star.chr.bam
samtools index ${OP}.star.chr.bam
rm ${OP}.header 


