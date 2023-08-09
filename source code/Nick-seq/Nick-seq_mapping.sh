#!/bin/bash
# Nick-seq_mapping


cutadapt -j 24 -a AGATCGGAAGAG -A AGATCGGAAGAG -u 3 ${sample}_cutadapt_1.fastq.gz ${sample}_cutadapt_2.fastq.gz -o ${sample}_1.fastq.gz -p ${sample}_2.fastq.gz
bowtie2 -p 24 -x $BOWTIE2INDEX -1 ${sample}_cutadapt_1.fastq.gz -2 ${sample}_cutadapt_2.fastq.gz  -S ${sample}.${SZ}.sam 
# if NT
	samtools view -@ 24 -f 64 -bSh ${sample}.${SZ}.sam  > ${sample}.${SZ}.bam 
# if TdT
	samtools view -@ 24 -f 128 -bSh ${sample}.${SZ}.sam  > ${sample}.${SZ}.bam 

samtools sort -@ 24 -o ${sample}.${SZ}.sorted.bam ${sample}.${SZ}.bam 
rm ${sample}.${SZ}.sam  ${sample}.${SZ}.bam 
bedtools bamtobed -i ${sample}.${SZ}.sorted.bam |bedtools intersect -v -a - -b $BLACKLIST >  ${sample}.${SZ}.sorted.bed 
rpm=$(wc -l ${sample}.${SZ}.sorted.bed |cut -f 1 -d ' '|awk '{printf "%f\n",$1/1000000.0}')
rpm=$(awk '{printf "%f\n",1/("'$rpm'")}' /sibcb1/wuweilab1/wuwei2022/Wuwei/Pipeline/tst.txt)
bedtools genomecov -bg  -scale $rpm -i ${sample}.${SZ}.sorted.bed -g $CHROMSIZE |bedtools sort -i stdin>  ${sample}.${SZ}.bedgraph
bedGraphToBigWig ${sample}.${SZ}.bedgraph $CHROMSIZE ${sample}.${SZ}.bw
rm ${sample}.${SZ}.bedgraph
date

# Use genomeCoverageBed to generate tabular files in GEO. We download the result file directly from GEO for downstream analysis following the guide of "https://github.com/BoCao2019/Nick-seq"






