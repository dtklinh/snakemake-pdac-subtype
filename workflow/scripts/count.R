log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

require(data.table)
library(Rsubread)

input <- snakemake@input
vectorInput <- unlist(input)
fc <- featureCounts(files=c(vectorInput),annot.ext=snakemake@config[["annotation"]],isGTFAnnotationFile=TRUE,GTF.featureType="exon",GTF.attrType="gene_id",nthreads=snakemake@threads,isPairedEnd=snakemake@config[["counting"]][["isPairedEnd"]])
tsv <- cbind(fc[["annotation"]]["GeneID"], fc[["counts"]])

samples = read.table(snakemake@config[["samples"]], header=TRUE, stringsAsFactors = FALSE)
setnames(tsv, c("GeneID", samples[["sample"]]))

write.table(tsv, "counts/all.tsv", sep="\t", quote=F, row.names=F)
