library("vsn")
library("pheatmap")
library("DESeq2")
library("BiocParallel")
library("RColorBrewer")
library("ggplot2")
library("ggrepel")

args<-commandArgs(TRUE)
print(args)

volcano_plot <- function(data,condition1,condition2,threshold){
resdata<-as.data.frame(data)
resdata$gene<-rownames(data)
resdata<-resdata[complete.cases(resdata),]
resdata$Legend<- ifelse(resdata$padj < 0.05 & abs(resdata$log2FoldChange)>1,"padj < 0.05 & FC > 2",
ifelse(resdata$padj < 0.05,"padj < 0.05",
ifelse(abs(resdata$log2FoldChange)>1, "FC > 2",
"Nonsignificant"
)))
top_genes<-resdata[resdata$Legend=="padj < 0.05 & FC > 2",]
top_genes<-top_genes[order(top_genes$baseMean,abs(top_genes$log2FoldChange),decreasing=TRUE),]
top_genes<-head(top_genes,threshold)
resdata$label <- ifelse(resdata$gene %in% top_genes$gene,resdata$gene,"")
cols <- c("Nonsignificant" = "darkgrey", "padj < 0.05 & FC > 2" = "#e69f00", "padj < 0.05" = "#56b4e9", "FC > 2" = "#009e73")

# Make a basic ggplot2 object

vol <- ggplot(resdata, aes(x = log2FoldChange, y = -log10(padj), colour = Legend))

# inserting manual colors as per color pallette  with term "scale_colour_manual(values = cols)" below

vol +   
  ggtitle(label = "Volcano Plot") +
  geom_point(size = 2.5, alpha = .99, na.rm = T)  +
  geom_text_repel(data=(resdata ), aes(x=resdata$log2FoldChange, y=-log10(padj),label=label),colour="black",point.padding = unit(1.6, "lines"))+
  scale_color_manual(name="Differential Expression",values=cols)+
  theme_bw(base_size = 14) + 
  theme(legend.position = "right") + 
  xlab(bquote(log[2](~.(condition1) / ~.(condition2)))) + 
  ylab(expression(-log[10]("pvalue"))) + # Change Y-Axis label
  geom_hline(yintercept = -log10(0.05), colour="#990000", linetype="dashed") + geom_vline(xintercept = 1, colour="#990000", linetype="dashed") + geom_vline(xintercept = -1, colour="#990000", linetype="dashed") + 
  scale_y_continuous(trans = "log1p",breaks=c(1,10,100,1000))

}

count_filename<-args[1] #count table
coldata_filename<-args[2] #sample description
outfolder <- args[3]
projectname <- args[4]
pthreads <- as.numeric(args[5])
condition1 <- args[6]
condition2 <-args[7]
threshold <- as.numeric(args[8])
register(MulticoreParam(pthreads))
transformation <- args[9] #all, ntd, rsd, vsd
batch<-args[10]
blind<-args[11]
timecourse<-args[12]
treatment<-args[13]
#read in file
setwd(outfolder)
cts<- as.matrix(read.delim(count_filename,header=TRUE, stringsAsFactors=TRUE,row.names="gene_id"))
coldata<-(read.delim(coldata_filename,header=TRUE, stringsAsFactors=TRUE,row.names="sample"))
rownames(coldata)<-colnames(cts)
### Set up Experiment Design
if (timecourse == "" & treatment == ""){
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
} else if (timecourse == "" & treatment != "") {
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~ condition + treatment)
} else if (timecourse != "" & treatment =="") {
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~ condition + timecourse + condition:timecourse)
}

if (batch=="sva"){
  library("sva")
  dat  <- counts(dds, normalized = TRUE)
  idx  <- rowMeans(dat) > 1
  dat  <- dat[idx, ]
  mod  <- model.matrix(~ dex, colData(dds))
  mod0 <- model.matrix(~   1, colData(dds))
  svseq <- svaseq(dat, mod, mod0, n.sv = 2)
  ddssva <- dds
  ddssva$SV1 <- svseq$sv[,1]
  ddssva$SV2 <- svseq$sv[,2]
  design(ddssva) <- ~ SV1 + SV2 + dex
  dds<-ddsva
} 
if (timecourse==""){
dds <- DESeq(dds,betaPrior=TRUE,parallel=TRUE, BPPARAM=MulticoreParam(pthreads))
normalized_counts=counts(dds, normalized=TRUE)
res <- results(dds,parallel=TRUE, BPPARAM=MulticoreParam(pthreads),contrast=c("condition",condition1,condition2))
resOrdered <- res[order(res$pvalue),]
summary(res)
} else {
  dds <- DESeq(dds,test="LRT", reduced= ~ condition + timecourse, betaPrior=TRUE,parallel=TRUE, BPPARAM=MulticoreParam(pthreads))
  normalized_counts=counts(dds, normalized=TRUE)
  res <- results(dds,parallel=TRUE, BPPARAM=MulticoreParam(pthreads),contrast=c("condition",condition1,condition2))
  resOrdered <- res[order(res$pvalue),]
  summary(res)
}

if (batch=="ruv") {
  library("RUVSeq")
  set <- newSeqExpressionSet(counts(dds))
  idx  <- rowSums(counts(set) > 5) >= 2
  set  <- set[idx, ]
  set <- betweenLaneNormalization(set, which="upper")
  not.sig <- rownames(res)[which(res$pvalue > .1)]
  empirical <- rownames(set)[ rownames(set) %in% not.sig ]
  set <- RUVg(set, empirical, k=2)
  pData(set)
  ddsruv <- dds
  ddsruv$W1 <- set$W_1
  ddsruv$W2 <- set$W_2
  design(ddsruv) <- ~ W1 + W2 + dex
  ddsruv <- DESeq(ddsruv,betaPrior=TRUE,parallel=TRUE, BPPARAM=MulticoreParam(pthreads))
  normalized_counts=counts(ddsruv, normalized=TRUE)
  res <- results(ddsruv,parallel=TRUE, BPPARAM=MulticoreParam(pthreads),contrast=c("condition",condition1,condition2))
  resOrdered <- res[order(res$pvalue),]
  summary(res)  
}

TE_only<-resOrdered[grepl(":",rownames(resOrdered)),]
RefSeq_only<-resOrdered[!rownames(resOrdered) %in% rownames(TE_only),]

write.table(resOrdered,"DESeq2_all.txt",row.names=TRUE,col.names=TRUE,sep="\\t",quote=FALSE)
write.table(TE_only,"DESeq2_TE_only.txt",row.names=TRUE,col.names=TRUE,sep="\\t",quote=FALSE)
write.table(RefSeq_only,"DESeq2_RefSeq_only.txt",row.names=TRUE,col.names=TRUE,sep="\\t",quote=FALSE)



#Count graphs
experiment <- ifelse(timecourse=FALSE & treatment=FALSE,"condition",
                     ifelse(timecourse=FALSE & treatment=TRUE,c("condition","treatment"),
                            ifelse(timecourse=TRUE & treatment=FALSE, c("condition","timecourse"))))

pdf("count_graph_all.pdf")
plotCounts(dds, gene=which.min(res$padj), intgroup=experiment)
dev.off()
# pdf("count_graph_TE_only.pdf")
# plotCounts(dds, gene=which.min(TE_only$padj), intgroup="condition")
# dev.off()
# pdf("count_graph_RefSeq_only.pdf")
# plotCounts(dds, gene=which.min(RefSeq_only$padj), intgroup="condition")
# dev.off()
#MA graphs
pdf("MA_plot_all.pdf")
plotMA(res, ylim=c(-2,2))
dev.off()
pdf("MA_plot_TE_only.pdf")
plotMA(TE_only, ylim=c(-2,2))
dev.off()
pdf("MA_plot_RefSeq_only.pdf")
plotMA(RefSeq_only, ylim=c(-2,2))
dev.off()

resSig<-res[ which(res$padj < .05), ]
TESig<- TE_only[which(TE_only$padj < 0.05),]
RefSeqSig<- RefSeq_only[which(RefSeq_only$padj < 0.05),]
if (nrow(TESig) > 0) {
  pdf("volcano_TE_only.pdf")
  volcano_plot(TE_only,condition1,condition2,threshold)
  dev.off()
}

if (nrow(RefSeqSig) > 0) {
  pdf("volcano_RefSeq_only.pdf")
  volcano_plot(RefSeq_only,condition1,condition2,threshold)
  dev.off()
}
if (nrow(resSig) > 0) { #if any significant genes
  #Volcano plots
  pdf("volcano_all.pdf")
  volcano_plot(res,condition1,condition2,threshold)
  dev.off()

  if (treatment == FALSE and timecourse == FALSE){ #if simple experiment
  if (transformation=="all"){ #if transformation not selected
  #Transformation Plots
  ntd <- normTransform(dds)
  vsd <- vst(dds, blind=FALSE)
  rld <- rlog(dds, blind=FALSE)
  
  pdf("meanSd_plot_all_ntd.pdf")
  meanSdPlot(assay(ntd))
  dev.off()
  pdf("meanSd_plot_all_vsd.pdf")
  meanSdPlot(assay(vsd))
  dev.off()
  pdf("meanSd_plot_all_rld.pdf")
  meanSdPlot(assay(rld))
  dev.off()
  
  select <- order(rowMeans(counts(dds,normalized=TRUE)),
                  decreasing=TRUE)[1:threshold]
  df <- as.data.frame(colData(dds)[,experiment])
  
  pdf("top100_heatmap_all_ntd.pdf")
  pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
           cluster_cols=FALSE, annotation_col=coldata)
  dev.off()
  pdf("top100_heatmap_all_vsd.pdf")
  pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
           cluster_cols=FALSE, annotation_col=coldata)
  dev.off()
  pdf("top100_heatmap_all_rld.pdf")
  pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
           cluster_cols=FALSE, annotation_col=coldata)
  dev.off()
  
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  
  pdf("sample_sample_distance_all_vsd.pdf")
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
  dev.off()
  
  pdf("pca_all_vsd.pdf")
  plotPCA(vsd, intgroup=experiment)
  dev.off()
  
  pdf("count_outliers_all.pdf")
  par(mar=c(8,5,2,2))
  boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
  dev.off()
  
  pdf("dispersion_all.pdf")
  plotDispEsts(dds)
  dev.off()
  
  metadata(res)$alpha
  metadata(res)$filterThreshold
  pdf("independent_filtering_rejections_all")
  plot(metadata(res)$filterNumRej, 
       type="b", ylab="number of rejections",
       xlab="quantiles of filter")
  lines(metadata(res)$lo.fit, col="red")
  abline(v=metadata(res)$filterTheta)
  dev.off()
  
  resNoFilt <- results(dds, independentFiltering=FALSE)
  addmargins(table(filtering=(res$padj < .1),
                   noFiltering=(resNoFilt$padj < .1)))
  
} else if (transformation =="vsd"){
  norm <- vst(dds, blind=FALSE)
} else if (transformation=="ntd"){
  norm <- normTransform(dds)
} else if (transformation=="rtd"){
  norm <- rlog(dds, blind=FALSE)
}
library("genefilter")
topVarGenes <- head(order(rowVars(assay(norm)), decreasing = TRUE), threshold)
mat  <- assay(norm)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(norm)[, experiment])
pheatmap(mat, annotation_col = anno)
}
    
}
}
}




# par(mfrow=c(2,2),mar=c(2,2,1,1))
# ylim <- c(-2.5,2.5)
# resGA <- results(dds, lfcThreshold=.5, altHypothesis="greaterAbs")
# resLA <- results(dds, lfcThreshold=.5, altHypothesis="lessAbs")
# resG <- results(dds, lfcThreshold=.5, altHypothesis="greater")
# resL <- results(dds, lfcThreshold=.5, altHypothesis="less")
# drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
# plotMA(resGA, ylim=ylim); drawLines()
# plotMA(resLA, ylim=ylim); drawLines()
# plotMA(resG, ylim=ylim); drawLines()
# plotMA(resL, ylim=ylim); drawLines()

# TE_only<-resOrdered[grepl("\\|",rownames(resOrdered)),]
# plotCounts(dds, gene=which.min(TE_only$padj), intgroup="condition")
# plotMA(TE_only, ylim=c(-2,2))



"""
with open(outfilepath,'w') as outfile:
outfile.write(Rscript)