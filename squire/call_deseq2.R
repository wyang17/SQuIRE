source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library("DESeq2")
library("BiocParallel")
register(MulticoreParam(4))
args<-commandArgs(TRUE)
print(args)

count_filename<-args[1] #count table
coldata_filename<-args[2] #sample description
strand<-args[3]
outfile <- args[4]
print(args)

#read in file
#count_filename="E:/Dropbox/Miguel/squire_call/squire_miguel_counttable.txt"
#coldata_filename="E:/Dropbox/Miguel/squire_call/squire_miguel_coldata.txt"
cts<- as.matrix(read.delim(count_filename,header=TRUE, stringsAsFactors=FALSE,row.names="gene_id"))
coldata<-(read.delim(coldata_filename,header=TRUE, stringsAsFactors=FALSE,row.names="sample"))

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds$condition <- factor(dds$condition, levels = c("control","treated"))
dds <- DESeq(dds)
res <- results(dds)
resLFC <- lfcShrink(dds, coef=2)
resOrdered <- res[order(res$pvalue),]
summary(res)

plotCounts(dds, gene=which.min(res$padj), intgroup="condition")
plotMA(res, ylim=c(-2,2))

TE_only<-resOrdered[grepl("\\|",rownames(resOrdered)),]
plotMA(TE_only, ylim=c(-2,2))
