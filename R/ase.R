library(GenomicFeatures)
library(ggplot2)
library(scales)
library(qvalue)

# Read ASE table
args=commandArgs(trailingOnly=TRUE)
as=read.table(args[1], header=T)

# Subset by depth
as=as[as$depth>20,]

# Plot VAF distribution
png(paste0(args[1], ".vaf.png"), height=600, width=600)
p1=ggplot(data=as, aes(x=af))
p1=p1 + geom_freqpoly(aes(y=..density..), binwidth=0.05)
p1=p1 + xlab("Variant Allele Frequency") + ylab("Density")
p1=p1 + scale_y_continuous(labels=comma) + ggtitle("RNA AF-Distribution of het. SNPs")
p1=p1 + scale_x_continuous(labels=comma, limits=c(0,1)) + geom_vline(xintercept=0.2, linetype="dashed") + geom_vline(xintercept=0.8, linetype="dashed")
p1
dev.off()
print(warnings())

# Multiple testing correction
m = nrow(as)
threshold = 0.05
as$bonferoni=(as$pvalue<threshold/m)  # Bonferoni
as$padj=(p.adjust(as$pvalue, method = 'fdr')<threshold)  # p-adjust fdr
as$qval=(qvalue(as$pvalue)$qvalues<threshold)  # q-values

print("Fraction allele-specific")
print(paste0("Bonferoni: ", mean(as$bonferoni)))
print(paste0("P-adjust FDR: ", mean(as$padj)))
print(paste0("q-values: ", mean(as$padj)))
library(ggplot2)
library(scales)

# All chromosomes
chrs = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22")
as = as[as$chr %in% chrs,]
as$chr = factor(as$chr, levels=chrs)

png(paste0(args[1], ".genome.png"), height=1200, width=1200)
p1 = ggplot(data=as, aes(x=pos, y=af))
p1 = p1 + geom_point(pch=21, size=0.5)
p1 = p1 + xlab("Chromosome")
p1 = p1 + ylab("Allele Frequency of het. SNPs")
p1 = p1 + scale_x_continuous(labels=comma) + scale_y_continuous(labels=comma, limits=c(0,1))
p1 = p1 + facet_wrap(~ chr, scales="free_x")
p1 = p1 + theme(axis.text.x = element_text(angle=45, hjust=1))
p1
dev.off()
print(warnings())

