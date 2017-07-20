library(GenomicFeatures)

exonTable = function(id) {
	  db=makeTxDbFromUCSC(genome=id, tablename="ccdsGene")
	  ex=reduce(exons(db), ignore.strand=T)
	  ex=keepStandardChromosomes(ex)
	  ex=ex[width(ex)>1,]
	  df=data.frame(chr=seqnames(ex), start=start(ex), end=end(ex), type="exonic")
	  df$chr=substring(as.character(df$chr), 4)	  
	  gz = gzfile(paste0("exonic.", id, ".bed.gz"), "w")
	  write.table(df, gz, quote=F, row.names=F, col.names=F, sep="\t")
	  close(gz)
}

exonTable("hg19")
