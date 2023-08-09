#!/bin/bash
# END-seq


bowtie -p 24 -n 3 -l 50 -k 1 ${BowtieIndex} -q ${File}.fastq.gz -S ${File}.sam --chunkmbs 200
samtools view -bS ${File}.sam > ${File}.bam
samtools sort -@ 24 -o ${File}.sorted.bam ${File}.bam
rm ${File}.bam
rm ${File}.sam
bedtools bamtobed -i ${File}.sorted.bam|bedtools intersect -v -a - -b ${Blacklist} > ${File}.sorted.bed
rpm=$(wc -l ${File}.sorted.bed|cut -f 1 -d ' '|awk '{printf "%f\n",$1/1000000.0}')
rpm=$(awk '{printf "%f\n",1/("'$rpm'")}' /sibcb1/wuweilab1/wuwei2022/Wuwei/Pipeline/tst.txt)
bedtools genomecov -bg -strand "+" -scale $rpm -i ${File}.sorted.bed -g ${chromSize}|bedtools sort -i stdin> ${File}.pos.bedgraph
bedtools genomecov -bg -strand "-" -scale $rpm -i ${File}.sorted.bed -g ${chromSize}|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin > ${File}.neg.bedgraph
bedtools genomecov -bg -5 -strand "+" -scale $rpm -i ${File}.sorted.bed -g ${chromSize}|bedtools sort -i stdin> ${File}.pos5.bedgraph
bedtools genomecov -bg -5 -strand "-" -scale $rpm -i ${File}.sorted.bed -g ${chromSize}|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin > ${File}.neg5.bedgraph

bedGraphToBigWig ${File}.pos.bedgraph ${chromSize} ${File}.pos.bw
bedGraphToBigWig ${File}.neg.bedgraph ${chromSize} ${File}.neg.bw
bedGraphToBigWig ${File}.pos5.bedgraph ${chromSize} ${File}.pos5.bw
bedGraphToBigWig ${File}.neg5.bedgraph ${chromSize} ${File}.neg5.bw

rm ${File}.pos.bedgraph
rm ${File}.neg.bedgraph
rm ${File}.pos5.bedgraph
rm ${File}.neg5.bedgraph
