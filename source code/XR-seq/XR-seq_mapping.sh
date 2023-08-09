#!/bin/bash
# XR-seq

python3 -m pip install --user cutadapt
## Small RNA 3' Adapter
cutadapt -j 24 -m 20 -a TGGAATTCTCGGGTGCCAAGG -o $DPATH/${File}_cutadapt.fastq.gz $DPATH/$File.fastq.gz
bowtie2 -p 24 -x $Bowtie2Index -U $DPATH/${File}_cutadapt.fastq.gz -S $DPATH/${File}_cutadapt.$SZ.sam --no-unal
samtools view -bSh -q 20 -@ 24$DPATH/ ${File}_cutadapt.$SZ.sam > $DPATH/${File}_cutadapt.$SZ.bam
rm ${File}_cutadapt.$SZ.sam
bedtools bamtobed -i $DPATH/${File}_cutadapt.$SZ.bam > $DPATH/${File}_cutadapt.$SZ.bed
sort -u -k1,1 -k2,2n -k3,3n $DPATH/${File}_cutadapt.$SZ.bed > $DPATH/${File}_cutadapt.$SZ.rmdup.sorted.bed
rm $DPATH/${File}_cutadapt.$SZ.bed
awk '{ if ($3-$2 == 26) print $0}' $DPATH/${File}_cutadapt.$SZ.rmdup.sorted.bed > $DPATH/${File}_cutadapt.$SZ.rmdup.sorted_26bp.bed
bedtools getfasta -fi $fa -bed $DPATH/${File}_cutadapt.$SZ.rmdup.sorted_26bp.bed -fo $DPATH/${File}_cutadapt.$SZ.rmdup.sorted_26bp.fa -s

rpm=$(wc -l $DPATH/${File}_cutadapt.$SZ.rmdup.sorted.bed|cut -f 1 -d ' '|awk '{printf "%f\n",$1/1000000.0}')
rpm=$(awk '{printf "%f\n",1/("'$rpm'")}' /sibcb1/wuweilab1/wuwei2022/Wuwei/Pipeline/tst.txt)

bedtools genomecov -bg -strand "+" -scale $rpm -i $DPATH/${File}_cutadapt.$SZ.rmdup.sorted.bed -g ${chromSize}|bedtools sort -i stdin> ${File}.pos.bedgraph
bedtools genomecov -bg -strand "-" -scale $rpm -i $DPATH/${File}_cutadapt.$SZ.rmdup.sorted.bed -g ${chromSize}|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin > ${File}.neg.bedgraph

bedGraphToBigWig ${File}.pos.bedgraph ${chromSize} ${File}.pos.bw
bedGraphToBigWig ${File}.neg.bedgraph ${chromSize} ${File}.neg.bw

rm ${File}.pos.bedgraph ${File}.neg.bedgraph
