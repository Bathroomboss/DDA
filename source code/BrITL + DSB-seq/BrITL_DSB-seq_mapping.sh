#!/bin/bash
# BrILT + DSB-seq


bowtie2 -p 24 --end-to-end --sensitive --score-min L,-1.5,-0.3 -x $BOWTIE2INDEX -U ${sample}.fastq.gz -S ${sample}.${SZ}.sam 
samtools view -@ 24 -q 13 -F 1024 -bSh ${sample}.${SZ}.sam > ${sample}.${SZ}.bam
samtools sort -@ 24 -o ${sample}.${SZ}.sorted.bam ${sample}.${SZ}.bam
rm ${sample}.${SZ}.sam  ${sample}.${SZ}.bam
bedtools bamtobed -i ${sample}.${SZ}.sorted.bam|bedtools intersect -v -a - -b $BLACKLIST  > ${sample}.${SZ}.sorted.bed

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