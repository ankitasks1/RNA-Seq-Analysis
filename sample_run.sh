for vars in sample*.fastq.gz
do
	base=$(basename ${vars} .fastq.gz)
	echo ${base}
	cutadapt -a AGATCGGAAGAG -m 30 -o ${base}_trimmed.fq.gz ${base}.fastq.gz
	
	#Alignment by Bowitei2: SE
	bowtie2 -q --phred33 -x ./ref_av/bowtie2/planaria -U ${base}_trimmed.fq.gz -S ${base}.sam
	#Convert SAM to BAM
	samtools view -S -b ${base}.sam > ${base}.bam

	#Extract Uniquely mapped reads
	samtools view -b -q 10 ${base}.bam > ${base}_uniq.bam

	#Sort BAM
	samtools sort -o ${base}_uniq_sorted.bam ${base}_uniq.bam
	samtools index ${base}_uniq_sorted.bam

	#Alignment by STAR
	STAR --runMode alignReads --runThreadN 10 --alignEndsType EndToEnd --genomeDir ./ref_av/starindex/ --readFilesIn ${base}_trimmed.fq.gz --outFileNamePrefix ${base}_star_ --outFilterMismatchNmax 5 --outFilterMultimapNmax 1 --outBAMsortingThreadN 12 --outSAMtype BAM SortedByCoordinate --outSAMattributes All --outSAMstrandField intronMotif --readFilesCommand gunzip -c

done
