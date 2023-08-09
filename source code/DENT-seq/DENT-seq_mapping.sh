#!/bin/bash
# DENT-seq_mapping
# Note: we did not perform the single NT analysis and just use the parameters mentioned in the publication

sample=SRX8786118-NickSeq_hg19_highMapQ
Bowtie2Index=/sibcb1/wuweilab1/wuwei-lab/Reference/Human/hg19/Bowtie2IndexplusrDNA/genome
SZ=hg19
bowtie2 -p 24  -X 800 -x ${Bowtie2Index} -1 ${sample}_1.fastq.gz -2 ${sample}_2.fastq.gz -S ${sample}.${SZ}.sam --no-discordant 
samtools view -bSh -@ 24  $File.$SZ.sam > $File.$SZ.bam
rm $File.$SZ.sam
samtools sort -@ 24 -o $File.$SZ.sorted.bam $File.$SZ.bam
rm $File.$SZ.bam
java -Xms20g -Xmx20g -XX:ParallelGCThreads=4 -jar /sibcb/program/install/picard/picard.jar MarkDuplicates I=$File.$SZ.sorted.bam O=$File.$SZ.rmDup.bam M=$File.rmDup_metrics.txt  REMOVE_DUPLICATES=true
samtools sort -@ 24 -o $File.$SZ.rmDup.sorted.bam $File.$SZ.rmDup.bam
bedtools bamtobed -i $File.$SZ.rmDup.sorted.bam |bedtools intersect -v -a - -b $BLACKLIST > $File.$SZ.sorted.bed
rpm=$(wc -l $File.$SZ.sorted.bed|cut -f 1 -d ' '|awk '{printf "%f\n",$1/1000000.0}')
rpm=$(awk '{printf "%f\n",1/("'$rpm'")}' /sibcb1/wuweilab1/wuwei2022/Wuwei/Pipeline/tst.txt)
bedtools genomecov -bg -scale $rpm -i $File.$SZ.sorted.bed -g $CHROMSIZE|bedtools sort -i stdin> $File.$SZ.bedgraph
bedGraphToBigWig $File.$SZ.bedgraph $CHROMSIZE $File.$SZ.bw
rm $File.$SZ.bedgraph


