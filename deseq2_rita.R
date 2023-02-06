library(DESeq2)
library(EDASeq)
library(EnhancedVolcano)
library(data.table)
library(RColorBrewer)

#et the count table using samtools
samtools idxstats ctl_rep1_uniq_sorted.bam > ctl_rep1.count.txt
samtools idxstats ctl_rep2_uniq_sorted.bam > ctl_rep2.count.txt
samtools idxstats ctl_rep3_uniq_sorted.bam > ctl_rep3.count.txt
samtools idxstats sample1_rep1_1_uniq_sorted.bam > sample1_rep1.count.txt
samtools idxstats sample1_rep2_1_uniq_sorted.bam > sample1_rep2.count.txt
samtools idxstats sample1_rep3_1_uniq_sorted.bam > sample1_rep3.count.txt

setwd("/Users/ankitverma/Documents/Archivio2/tutorial/teacheron/rita")
ctl1 <- read.table("ctl_rep1.count.txt", header = F, stringsAsFactors = F)
ctl2 <- read.table("ctl_rep2.count.txt", header = F, stringsAsFactors = F)
ctl3 <- read.table("ctl_rep3.count.txt", header = F, stringsAsFactors = F)
sample1 <- read.table("sample1_rep1.count.txt", header = F, stringsAsFactors = F)
sample2 <- read.table("sample1_rep2.count.txt", header = F, stringsAsFactors = F)
sample3 <- read.table("sample1_rep3.count.txt", header = F, stringsAsFactors = F)

#Give a column name
colnames(ctl1) <- c("feature", "ctl1len", "ctl1", "ctl1unmapped")
colnames(ctl2) <- c("feature", "ctl2len", "ctl2", "ctl2unmapped")
colnames(ctl3) <- c("feature", "ctl3len", "ctl3", "ctl3unmapped")
colnames(sample1) <- c("feature", "sample1len", "sample1", "sample1unmapped")
colnames(sample2) <- c("feature", "sample2len", "sample2", "sample2unmapped")
colnames(sample3) <- c("feature", "sample3len", "sample3", "sample3unmapped")

#Create a count matrix
countdata <- merge(ctl1, ctl2, by = "feature")
countdata <- merge(countdata, ctl3, by = "feature")
countdata <- merge(countdata, sample1, by = "feature")
countdata <- merge(countdata, sample2, by = "feature")
countdata <- merge(countdata, sample3, by = "feature")
countdata <- countdata[-1,]
colnames(countdata)
head(countdata)
rownames(countdata) <- countdata$feature


#Get only required columns
countData = as.matrix(countdata[,c(3,6,9,12,15,18)])
head(countData)
coldata <- read.table("coldata.txt" , header = TRUE, stringsAsFactors = FALSE)
rownames(coldata) <- colnames(countData)
head(coldata)
coldata <- coldata[,c("condition","replicate")]
coldata$condition <- factor(coldata$condition)
coldata$replicate <- factor(coldata$replicate)
all(rownames(coldata) == colnames(countData)) #should print TRUE
dds <- DESeqDataSetFromMatrix(countData =countData, colData = coldata, design = ~ condition)
dds
featureData <- data.frame(gene=rownames(countData))
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
ddscounts <- counts(dds, normalized=FALSE)
head(ddscounts)
dim(ddscounts)
colSums(ddscounts)
#View filtered count matrix: View(counts(dds))
#Normalization is the part of DESeq command: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
#Normalized separately:
ddsNorm <- estimateSizeFactors(dds)
sizeFactors(ddsNorm)
#Export normalized counts
ddsNormcounts <- counts(ddsNorm, normalized=TRUE)
head(ddsNormcounts)
dim(ddsNormcounts)
write.table(ddsNormcounts, "ddsNormcounts.txt", sep="\t", quote=F, col.names=NA)
#PCA plot
plotPCA(as.matrix(ddsNormcounts), labels=T, col =  c("darkgreen","darkgreen","darkgreen","darkred","darkred","darkred"))
#save as PCA_ddsNormcounts.svg

#DEG
#Normalization is the part of DESeq command: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
dds <- DESeq(dds)
#See size factors estimated by DESeq function: sizeFactors(dds)
res <- results(dds, contrast=c("condition","irradiated","untreated"))
res
summary(res)
res_sorted = res[order(rownames(res)),]
head(res_sorted)
sum(res$padj < 0.05, na.rm=TRUE)
res_sorted = res[order(rownames(res)),]
head(res_sorted)
summary(res_sorted)
write.table(res_sorted, "deseq2_results_res_sorted.txt", sep="\t", quote = FALSE, append = FALSE)
sum(res_sorted$padj < 0.05, na.rm=TRUE)
res_sorted$threshold <- as.logical(res_sorted$padj < 0.05)
deseq2_results_res0.05 <- res_sorted[which(res_sorted$threshold == TRUE),]
head(deseq2_results_res0.05)
dim(deseq2_results_res0.05)
summary(deseq2_results_res0.05)
write.table(deseq2_results_res0.05, "deseq2_results_res0.05_sorted.txt", sep="\t", quote = FALSE, append = FALSE)



#MA plot
plotMA(dds)
#plotMA_dds.svg
#Volcano plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                labSize = 3)
#Volcanoplot.svg

#Tag chromosomes
#head(res_sorted)
res_Deseq2_features <- data.frame(rownames(res_sorted))
colnames(res_Deseq2_features) <- "id"
head(res_Deseq2_features)
#Chromosome positions
chr.pos =  fread("WIX1annotationTable.txt",header=T)
head(chr.pos)
chr.pos.1 = merge(res_Deseq2_features, chr.pos, by="id", all.x=TRUE)
head(chr.pos.1)
dim(chr.pos.1)
head(res_sorted)
res_sorted["id"] <- rownames(res_sorted)
res_sorted <- data.frame(res_sorted)
res_sorted_gene <- merge(res_sorted, chr.pos.1, by="id", all.x=TRUE)
head(res_sorted_gene)
tail(res_sorted_gene)
dim(res_sorted_gene)
write.table(res_sorted_gene, "res_sorted_gene.txt", sep="\t", quote = FALSE, append = FALSE, row.names = F, col.names = F)


#Tag chromosomes for Deregulated genes
#head(deseq2_results_res0.05)
deseq2_results_res0.05 <- data.frame(deseq2_results_res0.05)
deseq2_results_res0.05["id"] <- rownames(deseq2_results_res0.05)
deseq2_results_res0.05_gene <- merge(deseq2_results_res0.05, chr.pos.1, by="id", all.x=TRUE)

head(deseq2_results_res0.05_gene)
tail(deseq2_results_res0.05_gene)
dim(deseq2_results_res0.05_gene)
write.table(deseq2_results_res0.05_gene, "deseq2_results_res0.05_gene.txt", sep="\t", quote = FALSE, append = FALSE, row.names = F, col.names = F)

deseq2_results_res0.05_gene.high <- deseq2_results_res0.05_gene[which(deseq2_results_res0.05_gene$baseMean > 200),]
deseq2_results_res0.05_gene.high.up <- head(deseq2_results_res0.05_gene.high[order(-deseq2_results_res0.05_gene.high$log2FoldChange),], 10)
deseq2_results_res0.05_gene.high.down <- head(deseq2_results_res0.05_gene.high[order(deseq2_results_res0.05_gene.high$log2FoldChange),], 10)

deseq2_results_res0.05_gene.high.up.id <- data.frame(id=deseq2_results_res0.05_gene.high.up$id)
deseq2_results_res0.05_gene.high.down.id <- data.frame(id=deseq2_results_res0.05_gene.high.down$id)
deseq2_results_res0.05_gene.high.up.down <- rbind.data.frame(deseq2_results_res0.05_gene.high.up.id, deseq2_results_res0.05_gene.high.down.id)
#Heatmap
library("pheatmap")
#Transform and scale
z_TddsNormcounts= scale(t(ddsNormcounts), center = TRUE, scale = TRUE)
z_ddsNormcounts <- data.frame(t(z_TddsNormcounts))
z_ddsNormcounts["id"] <- rownames(z_ddsNormcounts)
z_ddsNormcounts_top <- merge(z_ddsNormcounts, deseq2_results_res0.05_gene.high.up.down, by="id", all.y=T)
rownames(z_ddsNormcounts_top) <- z_ddsNormcounts_top$id
z_ddsNormcounts_top <- z_ddsNormcounts_top[,-1]
breaksList = seq(-1, 1, by = 0.1)

pheatmap(z_ddsNormcounts_top,
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(length(breaksList)),
         breaks = breaksList,
         cluster_cols=T, 
         cluster_rows=T,
         fontsize = 12,
         clustering_distance_cols = "euclidean", 
         clustering_method = "ward.D", 
         border_color = "black",
         cutree_cols = 1)
#pheatmap_z_ddsNormcounts_top.svg

deseq2_results_res0.05_gene_allup <- deseq2_results_res0.05_gene[which(deseq2_results_res0.05_gene$log2FoldChange > 1),]
deseq2_results_res0.05_gene_alldown <- deseq2_results_res0.05_gene[which(deseq2_results_res0.05_gene$log2FoldChange < -1),]
write.table(deseq2_results_res0.05_gene_allup, "deseq2_results_res0.05_gene_allup.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = T)
write.table(deseq2_results_res0.05_gene_alldown, "deseq2_results_res0.05_gene_alldown.txt", sep = "\t", quote = F, append = F, row.names = F, col.names = T)




