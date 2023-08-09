#!/bin/bash
# INDUCE-seq_mapping

date
method="INDUCE_SEQ"
platform="NEXTSEQ550"
threads="24"
echo $File

bwa mem -t "$threads" -M -R '@RG\tID:INDUCE_SEQ\tPL:NEXTSEQ550\tPU:0\tLB:INDUCE_SEQ\tSM:"$File"' $BWA_Index $File.fastq.gz > $File.$SZ.sam
awk '$6 !~ /[0-9]S/{print}' $File.$SZ.sam | samtools view -Shb -q 30 -F 256 -F 1024 -@ "$threads" - |
samtools sort -@ "$threads" -m 4G - -o $File.$SZ.sorted.bam
rm $File.$SZ.sam
bedtools bamtobed -i ${File}.$SZ.sorted.bam|bedtools intersect -v -a - -b ${Blacklist} > ${File}.$SZ.sorted.bed
rpm=$(wc -l ${File}.$SZ.sorted.bed|cut -f 1 -d ' '|awk '{printf "%f\n",$1/1000000.0}')
rpm=$(awk '{printf "%f\n",1/("'$rpm'")}' /sibcb1/wuweilab1/wuwei2022/Wuwei/Pipeline/tst.txt)
bedtools genomecov -bg -strand "+" -scale $rpm -i ${File}.$SZ.sorted.bed -g ${chromSize}|bedtools sort -i stdin> ${File}.pos.bedgraph
bedtools genomecov -bg -strand "-" -scale $rpm -i ${File}.$SZ.sorted.bed -g ${chromSize}|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin > ${File}.neg.bedgraph
bedtools genomecov -bg -5 -strand "+" -scale $rpm -i ${File}.$SZ.sorted.bed -g ${chromSize}|bedtools sort -i stdin> ${File}.pos5.bedgraph
bedtools genomecov -bg -5 -strand "-" -scale $rpm -i ${File}.$SZ.sorted.bed -g ${chromSize}|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin > ${File}.neg5.bedgraph

bedGraphToBigWig ${File}.pos.bedgraph ${chromSize} ${File}.pos.bw
bedGraphToBigWig ${File}.neg.bedgraph ${chromSize} ${File}.neg.bw
bedGraphToBigWig ${File}.pos5.bedgraph ${chromSize} ${File}.pos5.bw
bedGraphToBigWig ${File}.neg5.bedgraph ${chromSize} ${File}.neg5.bw

rm ${File}.pos.bedgraph
rm ${File}.neg.bedgraph
rm ${File}.pos5.bedgraph
rm ${File}.neg5.bedgraph

date
