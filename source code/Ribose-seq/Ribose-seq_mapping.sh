#!/bin/bash
# Ribose-seq

cutadapt -j 0 -e 0.1 -O 3  -m 20 --quality-cutoff 25  -a AGATCGGAAGAG -o ${File}_cutadapt.fastq.gz ${File}.fastq.gz 
umi_tools extract --stdin=${File}_cutadapt.fastq.gz --stdout=$File.processed.fastq.gz --extract-method=regex --bc-pattern='(?P<umi_1>.{8})'  
rm ${File}_cutadapt.fastq.gz

## uniq
bowtie -p 24 -m 1 $BOWTIEINDEX -q $File.processed.fastq.gz -S $File.uniq.$SZ.sam --chunkmbs 200 --no-unal 
samtools view -@ 24 -bSh $File.uniq.$SZ.sam > $File.uniq.$SZ.bam
rm $File.uniq.$SZ.sam
samtools sort -@ 24 -o $File.uniq.$SZ.sorted.bam $File.uniq.$SZ.bam
rm $File.uniq.$SZ.bam
samtools index $File.uniq.$SZ.sorted.bam
umi_tools dedup -I $File.uniq.$SZ.sorted.bam -S $File.uniq.$SZ.dedup.sorted.bam --mapping-quality=30
bedtools bamtobed -i $File.uniq.$SZ.dedup.sorted.bam > $File.uniq.$SZ.sorted.bed
rpm=$(wc -l $File.uniq.$SZ.sorted.bed|cut -f 1 -d ' '|awk '{printf "%f\n",$1/1000000.0}')
rpm=$(awk '{printf "%f\n",1/("'$rpm'")}' /sibcb1/wuweilab1/wuwei2022/Wuwei/Pipeline/tst.txt)
bedtools genomecov -bg -strand "+" -scale $rpm -i $File.uniq.$SZ.sorted.bed -g $CHROMSIZE|bedtools sort -i stdin> $File.uniq.$SZ.pos.bedgraph
bedtools genomecov -bg -strand "-" -scale $rpm -i $File.uniq.$SZ.sorted.bed -g $CHROMSIZE|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin> $File.uniq.$SZ.neg.bedgraph
bedGraphToBigWig $File.uniq.$SZ.pos.bedgraph $CHROMSIZE $File.uniq.$SZ.pos.bw
bedGraphToBigWig $File.uniq.$SZ.neg.bedgraph $CHROMSIZE $File.uniq.$SZ.neg.bw
rm $File.uniq.$SZ.pos.bedgraph
rm $File.uniq.$SZ.neg.bedgraph

bedtools genomecov -bg -5 -strand "+" -scale $rpm -i $File.uniq.$SZ.sorted.bed -g $CHROMSIZE|bedtools sort -i stdin> $File.uniq.$SZ.pos5.bedgraph
bedtools genomecov -bg -5 -strand "-" -scale $rpm -i $File.uniq.$SZ.sorted.bed -g $CHROMSIZE|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin> $File.uniq.$SZ.neg5.bedgraph
bedGraphToBigWig $File.uniq.$SZ.pos5.bedgraph $CHROMSIZE $File.uniq.$SZ.pos5.bw
bedGraphToBigWig $File.uniq.$SZ.neg5.bedgraph $CHROMSIZE $File.uniq.$SZ.neg5.bw
rm $File.uniq.$SZ.pos5.bedgraph
rm $File.uniq.$SZ.neg5.bedgraph




## all
bowtie -p 24 --all $BOWTIEINDEX -q $File.processed.fastq.gz -S $File.all.$SZ.sam --chunkmbs 200 --no-unal 
samtools view -@ 24 -bSh $File.all.$SZ.sam > $File.all.$SZ.bam
rm $File.all.$SZ.sam
samtools sort -@ 24 -o $File.all.$SZ.sorted.bam $File.all.$SZ.bam
rm $File.all.$SZ.bam
samtools index $File.all.$SZ.sorted.bam
umi_tools dedup -I $File.all.$SZ.sorted.bam -S $File.all.$SZ.dedup.sorted.bam --mapping-quality=30
bedtools bamtobed -i $File.all.$SZ.dedup.sorted.bam > $File.all.$SZ.sorted.bed
rpm=$(wc -l $File.all.$SZ.sorted.bed|cut -f 1 -d ' '|awk '{printf "%f\n",$1/1000000.0}')
rpm=$(awk '{printf "%f\n",1/("'$rpm'")}' /sibcb1/wuweilab1/wuwei2022/Wuwei/Pipeline/tst.txt)
bedtools genomecov -bg -strand "+" -scale $rpm -i $File.all.$SZ.sorted.bed -g $CHROMSIZE|bedtools sort -i stdin> $File.all.$SZ.pos.bedgraph
bedtools genomecov -bg -strand "-" -scale $rpm -i $File.all.$SZ.sorted.bed -g $CHROMSIZE|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin> $File.all.$SZ.neg.bedgraph
bedGraphToBigWig $File.all.$SZ.pos.bedgraph $CHROMSIZE $File.all.$SZ.pos.bw
bedGraphToBigWig $File.all.$SZ.neg.bedgraph $CHROMSIZE $File.all.$SZ.neg.bw
rm $File.all.$SZ.pos.bedgraph
rm $File.all.$SZ.neg.bedgraph

bedtools genomecov -bg -5 -strand "+" -scale $rpm -i $File.all.$SZ.sorted.bed -g $CHROMSIZE|bedtools sort -i stdin> $File.all.$SZ.pos5.bedgraph
bedtools genomecov -bg -5 -strand "-" -scale $rpm -i $File.all.$SZ.sorted.bed -g $CHROMSIZE|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin> $File.all.$SZ.neg5.bedgraph
bedGraphToBigWig $File.all.$SZ.pos5.bedgraph $CHROMSIZE $File.all.$SZ.pos5.bw
bedGraphToBigWig $File.all.$SZ.neg5.bedgraph $CHROMSIZE $File.all.$SZ.neg5.bw
rm $File.all.$SZ.pos5.bedgraph
rm $File.all.$SZ.neg5.bedgraph

