library(DESeq2)
library(edgeR)
library(ggplot2)

# data files 
coldata <- read.table("rnaseq_design.txt", sep="\t", row.names = 1, header=TRUE)
countdata <- read.table("rnaseq_gene_counts.csv", sep="\t", row.names = 1,  header=TRUE)
countdata <- as.matrix(countdata)
rnaseq_annotation <- read.table("rnaseq_annotation.txt", sep="\t", row.names = 1, header=TRUE)

#filtering low counts, prior counts of three were added to negate effects of 0 counts on Log2 fold change
keep <- addPriorCount(countdata, lib.size=NULL, offset=NULL, prior.count=3)
keep <- rowSums(cpm(countdata)>1) >= 2
countdata <- countdata[keep,]

#making the dds object
dds <- DESeqDataSetFromMatrix(countdata, coldata, design = ~ condition)

#re-levelling to make control the reference level
dds$condition <- relevel(dds$condition, ref="control") 
dds$condition

#diff expression analysis
dds <- DESeq(dds)
res <- results(dds)
res

#shortlist
inx1 <- res$padj < 0.05 & res$log2FoldChange > 1
inx2 <- res$padj < 0.05 & res$log2FoldChange < -1
res1 <- res[inx1,]
res2 <- res[inx2,]
res <- res[inx1 | inx2, ]
res <- res[complete.cases(res), ]

# merge two data frames by ID
res1 <- as.data.frame(res1)
res1_col1<- res1[0]
str_split_fixed(res1_col1, ":", 2)
res1 <- merge(x = res1, y = rnaseq_annotation[c('gene_name')], by=0)

res2 <- as.data.frame(res2)
res2
res2 <- merge(x = res2, y = rnaseq_annotation[c('gene_name')], by=0)

write.csv(res1, file="DESeq_results_padj0.05_lfc>1_upreg.csv")
write.csv(res2, file="DESeq_results_padj0.05_lfc<-1_downreg.csv")

#showing the log2 fold changes over mean of normalised counts
plotMA(res, main="DESeq2", ylim=c(-2,2))

#data quality assessment 
rld <- rlog(dds, blind=FALSE) #transformed count matrix

#transpose of the transformed count matrix to give sample-sample distances
sampleDists <- dist(t(assay(rld)))

#Heatmap showing the Euclidean distances between the samples as calculated from the regularized log transformation.
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, rld$condition)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap::pheatmap(sampleDistMatrix,
                   clustering_distance_rows=sampleDists,
                   clustering_distance_cols=sampleDists,
                   col=colors)

#principle component analysis based on distance matrix
pca <- plotPCA(rld, intgroup="condition")

#editing PCA plot
data <- plotPCA(rld, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition, shape=condition)) +
  geom_point(size=3) +
  geom_text(aes(label="control_vs_treated"), size=2, nudge_x=0, nudge_y = -0.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()

#looking at outliers (cook's distance)
assays(dds)[["cooks"]]

par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2) #consistent cook's distances among the samples

#looking at dispersion 
#The dispersion is a parameter describing how much the variance deviates from the mean
plotDispEsts(dds, main=name) #typical dispersion. Final estimates tending to the fitted estimates

#filtration criteria
#plot to show filtration method
#The mean of normalized counts provides an independent statistic for filtering the tests.
#dds will remove genes with low mean normalised counts (i.e. those on the left of the graph)
#thus leaving most genes with a low p-value
plot(res$baseMean+1, -log10(res$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))

metadata(res)$alpha
metadata(res)$filterThreshold
plot(metadata(res)$filterNumRej,
     type="b", ylab="number of rejections",
     xlab="quantiles of filter", main="RNA seq")
lines(metadata(res)$lo.fit, col="red")
abline(v=metadata(res)$filterTheta) #shows where the cutoff lies. 
#good level of stringency

upreg.genes <- as.array(res1$Row.names) 
downreg.genes <- as.array(res2$gene_name)
comparison <- Reduce(intersect, list(gene.data, downreg.genes))

biocLite("clusterProfiler")
library(clusterProfiler)
## shows transcriptional programs from gene clustering
ego.up <- enrichGO(upreg.genes, 'org.Hs.eg.db', keyType = 'ENSEMBL', ont="BP", pvalueCutoff=0.01)
dotplot(ego.up, showCategory=30)

ego <- enrichGO(downreg.genes, 'org.Hs.eg.db', keyType = 'ENSEMBL', ont="BP", pvalueCutoff=0.01)
dotplot(ego, showCategory=30)
                      
