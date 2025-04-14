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

# 
<img src="https://github.com/ankitasks1/RNA-Seq-Analysis/blob/main/Rplot01.png" width="128"/>
