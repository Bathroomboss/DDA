#!/bin/bash
# click-code-seq mapping

bowtie2 -x ${BOWTIE2INDEX} -U ${File}.fastq.gz -S ${File}.$SZ.sam --end-to-end  --trim5 8 
samtools view -bSh -q 20 -@ 24 ${File}.$SZ.sam > ${File}.$SZ.bam
samtools sort -@ 24 -o ${File}.$SZ.sorted.bam ${File}.$SZ.bam
samtools rmdup -s ${File}.$SZ.sorted.bam ${File}.$SZ.sorted.rmdup.bam
bamToBed -i ${File}.$SZ.sorted.rmdup.bam > ${File}.$SZ.sorted.bed

rpm=$(wc -l $File.$SZ.sorted.bed|cut -f 1 -d ' '|awk '{printf "%f\n",$1/1000000.0}')
rpm=$(awk '{printf "%f\n",1/("'$rpm'")}' /sibcb1/wuweilab1/wuwei2022/Wuwei/Pipeline/tst.txt)

bedtools genomecov -bg -strand "+" -scale $rpm -i $File.$SZ.sorted.bed -g $CHROMSIZE|bedtools sort -i stdin> $File.$SZ.pos.bedgraph
bedtools genomecov -bg -strand "-" -scale $rpm -i $File.$SZ.sorted.bed -g $CHROMSIZE|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin> $File.$SZ.neg.bedgraph

bedtools genomecov -bg -5 -strand "+" -scale $rpm -i $File.$SZ.sorted.bed -g $CHROMSIZE|bedtools sort -i stdin> $File.$SZ.pos5.bedgraph
bedtools genomecov -bg -5 -strand "-" -scale $rpm -i $File.$SZ.sorted.bed -g $CHROMSIZE|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin> $File.$SZ.neg5.bedgraph

bedGraphToBigWig $File.$SZ.pos.bedgraph $CHROMSIZE $File.$SZ.pos.bw
bedGraphToBigWig $File.$SZ.neg.bedgraph $CHROMSIZE $File.$SZ.neg.bw
bedGraphToBigWig $File.$SZ.pos5.bedgraph $CHROMSIZE $File.$SZ.pos5.bw
bedGraphToBigWig $File.$SZ.neg5.bedgraph $CHROMSIZE $File.$SZ.neg5.bw

rm *.bedgraph