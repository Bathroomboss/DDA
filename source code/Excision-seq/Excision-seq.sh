#!/bin/bash

# not only sacCer3 but also E.coli
# For Predigestion

bowtie -p 24 -m 1 $BOWTIEINDEX -q $File.fastq.gz -S $File.uniq.$SZ.sam --chunkmbs 200 --no-unal
samtools view -@ 24 -bSh $File.uniq.$SZ.sam > $File.uniq.$SZ.bam
rm $File.uniq.$SZ.sam
samtools sort -@ 24 -o $File.uniq.$SZ.sorted.bam $File.uniq.$SZ.bam
rm $File.uniq.$SZ.bam

bedtools bamtobed -i $File.uniq.$SZ.sorted.bam > $File.uniq.$SZ.sorted.bed
rpm=$(wc -l $File.uniq.$SZ.sorted.bed|cut -f 1 -d ' '|awk '{printf "%f\n",$1/1000000.0}')
rpm=$(awk '{printf "%f\n",1/("'$rpm'")}' /sibcb1/wuweilab1/wuwei2022/Wuwei/Pipeline/tst.txt)
bedtools genomecov -bg -scale $rpm -i $File.uniq.$SZ.sorted.bed -g $CHROMSIZE|bedtools sort -i stdin> $File.$SZ.bedgraph
bedtools genomecov -bg -5 -scale $rpm -i $File.uniq.$SZ.sorted.bed -g $CHROMSIZE|bedtools sort -i stdin> $File.$SZ.5.bedgraph
bedGraphToBigWig $File.$SZ.bedgraph $CHROMSIZE $File.$SZ.bw
bedGraphToBigWig $File.$SZ.5.bedgraph $CHROMSIZE $File.$SZ.5.bw
rm $File.$SZ.bedgraph $File.$SZ.5.bedgraph

# For Postdigestion

bowtie -p 24 -m 1 $BOWTIEINDEX -q $File.fastq.gz -S $File.uniq.$SZ.sam --chunkmbs 200 --no-unal
samtools view -@ 24 -bSh $File.uniq.$SZ.sam > $File.uniq.$SZ.bam
rm $File.uniq.$SZ.sam
samtools sort -@ 24 -o $File.uniq.$SZ.sorted.bam $File.uniq.$SZ.bam
rm $File.uniq.$SZ.bam

bedtools bamtobed -i $File.uniq.$SZ.sorted.bam > $File.uniq.$SZ.sorted.bed
rpm=$(wc -l $File.uniq.$SZ.sorted.bed|cut -f 1 -d ' '|awk '{printf "%f\n",$1/1000000.0}')
rpm=$(awk '{printf "%f\n",1/("'$rpm'")}' /sibcb1/wuweilab1/wuwei2022/Wuwei/Pipeline/tst.txt)
bedtools genomecov -bg -scale $rpm -i $File.uniq.$SZ.sorted.bed -g $CHROMSIZE|bedtools sort -i stdin> $File.$SZ.bedgraph
bedGraphToBigWig $File.$SZ.bedgraph $CHROMSIZE $File.$SZ.bw
rm $File.$SZ.bedgraph

