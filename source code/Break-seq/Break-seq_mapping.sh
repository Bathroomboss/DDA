#!/bin/bash
# Break-seq


bowtie2 -p 24 -x $BT2_index -1 ${File}_1.fastq.gz -2 ${File}_2.fastq.gz -S ${File}.$SZ.sam --no-unal
samtools view -bSh -@ 24 ${File}.$SZ.sam > $File.$SZ.bam
samtools sort -@ 24 -o $File.$SZ.sorted.bam $File.$SZ.bam
rm ${File}.$SZ.sam $File.$SZ.bam
java -Xms20g -Xmx20g -XX:ParallelGCThreads=4 -jar /sibcb/program/install/picard/picard.jar MarkDuplicates I=$File.$SZ.sorted.bam O=$File.$SZ.rmDup.bam M=$File.rmDup_metrics.txt REMOVE_DUPLICATES=true
samtools sort -@ 24 -o $File.$SZ.sorted.rmDup.bam $File.$SZ.rmDup.bam
bedtools bamtobed -i $File.$SZ.sorted.rmDup.bam |bedtools intersect -v -a - -b $BLACKLIST > $File.$SZ.sorted.bed
rpm=$(wc -l $File.$SZ.sorted.bed|cut -f 1 -d ' '|awk '{printf "%f\n",$1/1000000.0}')
rpm=$(awk '{printf "%f\n",1/("'$rpm'")}' /sibcb1/wuweilab1/wuwei2022/Wuwei/Pipeline/tst.txt)
bedtools genomecov -bg -scale $rpm -i $File.$SZ.sorted.bed -g $CHROMSIZE|bedtools sort -i stdin> $File.$SZ.bedgraph
bedGraphToBigWig $File.$SZ.bedgraph $CHROMSIZE $File.$SZ.bt2.bw
rm $File.$SZ.bedgraph
date

