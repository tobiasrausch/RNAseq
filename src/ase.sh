#!/bin/bash

if [ $# -ne 4 ]
then
    echo ""
    echo "Usage: $0 <genome.fa> <aligned.bam> <input.vcf.gz> <outprefix>"
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
HG=${1}
BAM=${2}
VCF=${3}
OP=${4}

# Find sample name
BAMID="NA"
if [ `samtools view -H ${BAM} | grep -P "\tSM:" | wc -l` -eq 1 ]
then
    BAMID=`samtools view -H ${BAM} | grep -P "\tSM:" | sed 's/^.*SM://' | sed 's/\t.*$//'`
else
    BAMID=`echo ${BAM} | sed 's/^.*\///' | sed 's/\..*$//'`
fi
SIZE=`echo ${BAMID} | awk '{print length($1);}'`

# Find matching VCF sample
VCFID="NA"
for I in `seq ${SIZE} -1 1`
do
    PREFIX=`echo ${BAMID} | cut -c 1-${I}`
    if [ `bcftools view ${VCF} | grep -m 1 "^#CHROM" | cut -f 10- | tr '\t' '\n' | grep "^${PREFIX}" | wc -l` -eq 1 ]
    then
	VCFID=`bcftools view ${VCF} | grep -m 1 "^#CHROM" | cut -f 10- | tr '\t' '\n' | grep "^${PREFIX}"`
	break
    fi
done

# Matching samples found?
if [ ${VCFID} == "NA" ]
then
    echo "No matching sample IDs found!"
    echo ${BAMID} ${VCFID}
    exit;
else
    echo "Using ${BAMID} and ${VCFID}."
fi

# Subset to sample of interest and bi-allelic only
bcftools view -s ${VCFID} -m2 -M2 --types snps,indels ${VCF} | sed 's/^##contig=<ID=chr/##contig=<ID=/' | sed 's/^chr//' | bcftools annotate -O b -o ${OP}.asein.bcf -x INFO,^FORMAT/GT -
bcftools index ${OP}.asein.bcf

# Haplotag BAM file
alfred split -i -r ${HG} -v ${OP}.asein.bcf -p ${OP}.haplotagged.bam -s ${VCFID} ${BAM}

# Create ASE table
alfred ase -p -f -r ${HG} -v ${OP}.asein.bcf -s ${VCFID} -a ${OP}.ase.tsv.gz ${OP}.haplotagged.bam

# Annotate ASE table
zcat ${OP}.ase.tsv.gz | tail -n +2 | awk '{print $1"\t"$2"\t"($2+1)"\tID"NR;}' > ${OP}.ase.bed
alfred annotate -g ${BASEDIR}/../gtf/Homo_sapiens.GRCh37.75.gtf.gz -o ${OP}.anno.bed ${OP}.ase.bed
rm ${OP}.ase.bed
paste <(zcat ${OP}.ase.tsv.gz | head -n 1) <(head -n 1 ${OP}.anno.bed | cut -f 5,6) > ${OP}.ase.tsv
paste <(zcat ${OP}.ase.tsv.gz | tail -n +2 | sed 's/\tNA\tNA\tNA/\tNA\tNA/' | sort -k1,1V -k2,2n) <(tail -n +2 ${OP}.anno.bed | sort -k1,1V -k2,2n | cut -f 5,6) >> ${OP}.ase.tsv
rm ${OP}.anno.bed ${OP}.ase.tsv.gz
gzip ${OP}.ase.tsv

# Clean-up
rm ${OP}.asein.bcf ${OP}.asein.bcf.csi ${OP}.haplotagged.bam ${OP}.haplotagged.bam.bai

# Deactivate environment
source deactivate
