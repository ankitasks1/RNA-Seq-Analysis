# sra tool kit
fastq-dump --split-files SRR13510468 --gzip

# QC
fastqc SRR13510468.fastq.gz

# Trimming bad quality/adapter containing  reads
java -jar trimmomatic-0.39.jar PE input_R1.fq.gz input_R2.fq.gz R1_paired.fastq.gz  R1_unpaired.fastq.gz  R2_paired.fastq.gz  R2_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36

# QC
fastqc trimmed.fastq.gz

# Generate reference genome Index
STAR-2.5.4a/bin/Linux_x86_64/STAR --runThreadN 10 --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles mm10.fa --sjdbGTFfile gencode.vM20.chr_patch_hapl_scaff.annotation.gtf --sjdbOverhang 

# Alignment
STAR-2.5.4a/bin/Linux_x86_64_static/STAR --runMode alignReads --runThreadN 10 --genomeDir /home/ankitv/ref_av/gencodes/gencode_M20/ --sjdbGTFfile gencode.vM20.chr_patch_hapl_scaff.annotation.gtf --sjdbOverhang 124 --readFilesIn R1_paired.fastq.gz R2_paired.fastq.gz --outFileNamePrefix sample1_ --outFilterMismatchNmax 3 --outFilterMultimapNmax 1 --outBAMsortingThreadN 12 --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outSAMattributes All --outSAMstrandField intronMotif --readFilesCommand gunzip -c

# Index bam
samtools index *.sortedByCoord.out.bam

# Bedgraph for visualization
STAR-2.5.4a/bin/Linux_x86_64_static/STAR --runMode inputAlignmentsFromBAM --runThreadN 10 --inputBAMfile sample1_Aligned.sortedByCoord.out.bam --outWigType bedGraph --outWigNorm RPM --outWigStrand Stranded --outFileNamePrefix sample1_

# Featurecounts
featureCounts -t exon -g gene_id -s 2 -a /export/home/ankitv/ref_av/gencodes/gencode_ref_star/mouse.gencode.vM1.annotation.gtf -T 12 -o Sample_1_star-featureCounts.txt ./../Sample_1_Aligned.sortedByCoord.out.bam


# DESeq2
library(tidyverse)

library(readxl)

library(DESeq2)

library(edgeR)

library(clusterProfiler)

library(org.At.tair.db)

library(ggpubr)

setwd("~/Downloads/kalki")

#### read the count matrix

ctx <- read.table("~/Downloads/kalki/appliedgenomics_ankit/AppliedG_HW3/AppliedG_HW3/TARGET_counts.txt", header = T)

dim(ctx)

head(ctx)

rownames(ctx) <- gsub("gene:","",ctx$Geneid)

ctx <- ctx[,-1]

colnames(ctx) <- gsub("_Aligned.sortedByCoord.out.bam","", colnames(ctx))

colnames(ctx) <- gsub("..HWTKNDRX3_n01_","", colnames(ctx))

ctx <- ctx[,-c(17:18)]

#### read the metadata

coldata <- read.table("~/Downloads/kalki/coldata.txt", header = T)

coldata$Sample <- gsub("_Aligned.sortedByCoord.out.bam","", coldata$Sample)

coldata$Sample <- gsub("..HWTKNDRX3_n01_","", coldata$Sample)

rownames(coldata) <- coldata$Sample

coldata$TF <- factor(coldata$TF)

coldata$Replicate <- factor(coldata$Replicate)

all(rownames(coldata) == colnames(ctx))

### deseq object

dds <- DESeqDataSetFromMatrix(countData = ctx, colData = coldata, design = ~TF)

dds

keep <- rowSums(counts(dds) >= 5) >=3

dds <- dds[keep,]

ddsNorm <- estimateSizeFactors(dds)

sizeFactors(ddsNorm)

vst <- vst(dds, blind=TRUE)

pca_plot <-plotPCA(vst, intgroup=c("TF","Replicate"))

pca_plot

### DE analysis with DESeq2

dds <- DESeq2::DESeq(dds)


### de_temp are differentially expressed genes

### temp all genes 

<code>

  diff_deseq2_sig_list <- list()
fc = 1
for (i in c("NFYA6", "NLP5", "PDF2", "NFYA5")){
  print(i)
  temp <- results(dds, contrast=c("TF",i,"EV"))
  temp$threshold <- as.logical(temp$padj < 0.05)
  de_temp <- temp[which(temp$threshold == TRUE),]
  de_temp <- data.frame(de_temp)
  de_temp_sig <- de_temp %>% dplyr::filter((log2FoldChange > fc) | (log2FoldChange < -fc))
  de_temp_up <- de_temp %>% dplyr::filter((log2FoldChange > fc))
  de_temp_down <- de_temp %>% dplyr::filter((log2FoldChange < -fc))
  diff_deseq2_sig_list[[paste0(i, "_", "EV")]][["all"]] <- temp
  diff_deseq2_sig_list[[paste0(i, "_", "EV")]][["de"]] <- de_temp
  diff_deseq2_sig_list[[paste0(i, "_", "EV")]][["de_sig_all"]] <- de_temp_sig
  diff_deseq2_sig_list[[paste0(i, "_", "EV")]][["de_sig_up"]] <- de_temp_up
  diff_deseq2_sig_list[[paste0(i, "_", "EV")]][["de_sig_down"]] <- de_temp_down
}


### MA plot
deseq2_MA_plot_list <- list()
for (contrasts in names(diff_deseq2_sig_list)){
  print(contrasts)
  fdr = 0.05
  fc = 1
  deseq2_MA_plot_list[[contrasts]] <- ggmaplot(diff_deseq2_sig_list[[contrasts]]$all, fdr = fdr, fc = fc, size = 0.3, palette = c("#D22B2B", "#1465AC", "darkgray"),
                                                 genenames = rownames(diff_deseq2_sig_list[[contrasts]]$all), legend = "top", top = 20, font.label = c("bold", 5),
                                                 font.legend = "bold", font.main = "bold", ggtheme = ggplot2::theme_minimal())
  
}
</code>
##### ggmaplot is part of ggpubr package

### clusterprofiler analysis
<code>
cp_analysis_list <- list()
for (clist in names(diff_deseq2_sig_list)){
  print(clist)
  gene_list <- diff_deseq2_sig_list[[clist]]$de_sig_all$log2FoldChange
  names(gene_list) <- rownames(diff_deseq2_sig_list[[clist]]$de_sig_all)
  gene_list = sort(gene_list, decreasing = TRUE)
  message("running GO analysis")
  #### GO analysis
  gse <- gseGO(geneList=gene_list, 
               ont ="ALL", 
               keyType = "TAIR", 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = org.At.tair.db, 
               pAdjustMethod = "BH")
  message("running dotplot")
  gse_dotplot <- dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
  
  message("running ridge plot")
  gse_ridgeplot <- ridgeplot(gse) + labs(x = "enrichment distribution")
  
  message("running GSEA plot")
  gse_gseaplot <- gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)
  
  message("storing info and plots")
  cp_analysis_list[[clist]][["gse"]] <- gse
  cp_analysis_list[[clist]][["dotplot"]] <- gse_dotplot
  cp_analysis_list[[clist]][["ridgeplot"]] <- gse_ridgeplot
  cp_analysis_list[[clist]][["gseaplot"]] <- gse_gseaplot
}
</code>
# 
<img src="https://github.com/ankitasks1/RNA-Seq-Analysis/blob/main/Rplot01.png" width="128"/>
