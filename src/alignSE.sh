#!/bin/bash

if [ $# -ne 3 ]
then
    echo ""
    echo "Usage: $0 <read.fq.gz> <output prefix> <strand [0|1|2]>"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

# Activate environment
export PATH=${BASEDIR}/../bin/bin:${PATH}
source activate ${BASEDIR}/../bin/envs/rna


##############
# Remove with next alfred release
export PATH=/g/solexa/home/rausch/scripts/cpp/alfred/bin:${PATH}
#############


# CMD params
THREADS=4
FQ=${1}
OP=${2}
STRAND=${3}
HG=${BASEDIR}/../genome/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa

# Fastqc
mkdir -p ${OP}.read.fastqc/
fastqc -t ${THREADS} -o ${OP}.read.fastqc/ ${FQ}

# STAR alignment
STAR --runThreadN ${THREADS} --outFileNamePrefix ${OP}.star --outTmpDir ${OP}.tmpSTAR --genomeDir ${BASEDIR}/../genome/ --readFilesIn ${FQ} --readFilesCommand zcat
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
alfred count_rna -s ${STRAND} -g ${BASEDIR}/../gtf/Homo_sapiens.GRCh37.75.gtf.gz -o ${OP}.gene.count ${OP}.star.bam
alfred count_rna -s ${STRAND} -n fpkm -g ${BASEDIR}/../gtf/Homo_sapiens.GRCh37.75.gtf.gz -o ${OP}.gene.fpkm ${OP}.star.bam
alfred count_rna -s ${STRAND} -n fpkm_uq -g ${BASEDIR}/../gtf/Homo_sapiens.GRCh37.75.gtf.gz -o ${OP}.gene.fpkm_uq ${OP}.star.bam

# Create browser tracks
alfred tracks -o ${OP}.bedGraph.gz ${OP}.star.bam
igvtools totdf ${OP}.bedGraph.gz ${OP}.tdf hg19

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

# Deactivate environment
source deactivate
