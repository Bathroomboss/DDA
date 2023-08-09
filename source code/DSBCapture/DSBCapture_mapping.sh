#!/bin/bash
# DSBCapture

if [ ${END} == "SE" ] ;then
        cutadapt -j 24 -O 3 -a AGATCGGAAGAGC -o ${sample}_cutadapt.fastq.gz ${sample}.fastq.gz
else
        cutadapt -j 24 -O 3 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -o ${sample}_cutadapt_1.fastq.gz -p ${sample}_cutadapt_2.fastq.gz ${sample}_1.fastq.gz ${sample}_2.fastq.gz
fi

# 2 mapping
if [ ${END} == "SE" ] ;then
        bwa mem -M -t 24 $BWA_Index ${sample}_cutadapt.fastq.gz > ${sample}.${SZ}.sam
else
        bwa mem -M -t 24 $BWA_Index ${sample}_cutadapt_1.fastq.gz ${sample}_cutadapt_2.fastq.gz > ${sample}.${SZ}.sam
fi
samtools view -@ 24 -F 2820 -bSh ${sample}.${SZ}.sam  > ${sample}.${SZ}.bam
rm ${sample}.${SZ}.sam
samtools sort -@ 24 -o ${sample}.${SZ}.sorted.bam ${sample}.${SZ}.bam
rm ${sample}.${SZ}.bam
if [ ${END} == "SE" ] ;then
        bedtools bamtobed -i ${sample}.${SZ}.sorted.bam|bedtools intersect -v -a - -b $BLACKLIST  > ${sample}.${SZ}.sorted.bed
else
        java -Xms20g -Xmx20g -XX:ParallelGCThreads=8 -jar /sibcb/program/install/picard/picard.jar MarkDuplicates I=${sample}.${SZ}.sorted.bam O=${sample}.${SZ}.rmDup.sorted.bam M=${sample}.rmDup_metrics.txt ASO=coordinate REMOVE_DUPLICATES=true
        samtools view -f 65 -@ 24 -q 10 -F 1024 -bSh ${sample}.${SZ}.rmDup.sorted.bam >  ${sample}.${SZ}.rmDup.sorted.R1.bam
        bedtools bamtobed -i ${sample}.${SZ}.rmDup.sorted.R1.bam|bedtools intersect -v -a - -b $BLACKLIST  > ${sample}.${SZ}.sorted.bed
fi
rpm=$(wc -l ${sample}.${SZ}.sorted.bed|cut -f 1 -d ' '|awk '{printf "%f\n",$1/1000000.0}')
rpm=$(awk '{printf "%f\n",1/("'$rpm'")}' /sibcb1/wuweilab1/wuwei2022/Wuwei/Pipeline/tst.txt)
bedtools genomecov -bg -strand "+" -scale $rpm -i ${sample}.${SZ}.sorted.bed -g $CHROMSIZE|bedtools sort -i stdin> ${sample}.${SZ}.pos.bedgraph
bedtools genomecov -bg -strand "-" -scale $rpm -i ${sample}.${SZ}.sorted.bed -g $CHROMSIZE|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin > ${sample}.${SZ}.neg.bedgraph

bedtools genomecov -5 -bg -strand "+" -scale $rpm -i ${sample}.${SZ}.sorted.bed -g $CHROMSIZE|bedtools sort -i stdin> ${sample}.${SZ}.pos5.bedgraph
bedtools genomecov -5 -bg -strand "-" -scale $rpm -i ${sample}.${SZ}.sorted.bed -g $CHROMSIZE|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin > ${sample}.${SZ}.neg5.bedgraph

bedGraphToBigWig ${sample}.${SZ}.pos.bedgraph $CHROMSIZE ${sample}.${SZ}.pos.bw
bedGraphToBigWig ${sample}.${SZ}.pos5.bedgraph $CHROMSIZE ${sample}.${SZ}.pos5.bw
bedGraphToBigWig ${sample}.${SZ}.neg.bedgraph $CHROMSIZE ${sample}.${SZ}.neg.bw
bedGraphToBigWig ${sample}.${SZ}.neg5.bedgraph $CHROMSIZE ${sample}.${SZ}.neg5.bw

rm ${sample}.${SZ}.pos.bedgraph ${sample}.${SZ}.pos5.bedgraph ${sample}.${SZ}.neg.bedgraph ${sample}.${SZ}.neg5.bedgraph

bedtools genomecov -bg -scale $rpm -i ${sample}.${SZ}.sorted.bed -g $CHROMSIZE|bedtools sort -i stdin> ${sample}.${SZ}.bedgraph
bedGraphToBigWig ${sample}.${SZ}.bedgraph $CHROMSIZE ${sample}.${SZ}.bw
rm ${sample}.${SZ}.bedgraph
date



