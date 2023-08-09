#!/bin/bash

if [ ${END} == "SE" ];then
        bowtie -p 24 -n 2 -l 50 -m 1 $BOWTIEINDEX -q $DPATH/${File}.fastq.gz -S $DPATH/$File.$SZ.sam --chunkmbs 200 --no-unal
else
        bowtie -p 24 -n 2 -l 50 -m 1 $BOWTIEINDEX -1 $DPATH/${File}_1.fastq.gz -2 $DPATH/${File}_2.fastq.gz -S $DPATH/$File.$SZ.sam --chunkmbs 200 --no-unal
fi

samtools view -bSh -@ 24  $DPATH/$File.$SZ.sam > $DPATH/$File.$SZ.bam
rm $DPATH/$File.$SZ.sam
samtools sort -@ 24 -o $DPATH/$File.$SZ.sorted.bam $DPATH/$File.$SZ.bam
rm $DPATH/$File.$SZ.bam
java -Xms20g -Xmx20g -XX:ParallelGCThreads=4 -jar /sibcb/program/install/picard/picard.jar MarkDuplicates I=$DPATH/$File.$SZ.sorted.bam O=$DPATH/$File.$SZ.rmDup.bam M=$DPATH/$File.rmDup_metrics.txt REMOVE_DUPLICATES=true
samtools sort -@ 24 -o $DPATH/$File.$SZ.rmDup.sorted.bam $DPATH/$File.$SZ.rmDup.bam
bedtools bamtobed -i $DPATH/$File.$SZ.rmDup.sorted.bam |bedtools intersect -v -a - -b $BLACKLIST > $DPATH/$File.$SZ.sorted.bed

rm $DPATH/$File.$SZ.sorted.bam $DPATH/$File.$SZ.rmDup.bam

rpm=$(wc -l $DPATH/$File.$SZ.sorted.bed|cut -f 1 -d ' '|awk '{printf "%f\n",$1/1000000.0}')
rpm=$(awk '{printf "%f\n",1/("'$rpm'")}' /sibcb1/wuweilab1/wuwei2022/Wuwei/Pipeline/tst.txt)

if [ $STRAND == "FALSE" ];then
        bedtools genomecov -bg -scale $rpm -i $DPATH/$File.$SZ.sorted.bed -g $CHROMSIZE|bedtools sort -i stdin> $DPATH/$File.$SZ.bedgraph
        bedGraphToBigWig $DPATH/$File.$SZ.bedgraph $CHROMSIZE $DPATH/$File.$SZ.bw
        rm $DPATH/$File.$SZ.bedgraph
else
        bedtools genomecov -bg -strand "+" -scale $rpm -i $DPATH/$File.$SZ.sorted.bed -g $CHROMSIZE|bedtools sort -i stdin> $DPATH/$File.$SZ.pos.bedgraph
        bedtools genomecov -bg -strand "-" -scale $rpm -i $DPATH/$File.$SZ.sorted.bed -g $CHROMSIZE|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin> $DPATH/$File.$SZ.neg.bedgraph
        bedGraphToBigWig $DPATH/$File.$SZ.pos.bedgraph $CHROMSIZE $DPATH/$File.$SZ.pos.bw
        bedGraphToBigWig $DPATH/$File.$SZ.neg.bedgraph $CHROMSIZE $DPATH/$File.$SZ.neg.bw
        rm $DPATH/$File.$SZ.pos.bedgraph
        rm $DPATH/$File.$SZ.neg.bedgraph
fi

### NOTE!!!!!NOTE!!!NOTE!!!NOTE!!!NOTE!!!   Following code are prepared for PAIRED and strand specific ChIP-seq (such as TCR-seq, SSDS) to produce bigwig file.  NOTE!!!!!NOTE!!!NOTE!!!NOTE!!!NOTE!!!  

source /sibcb/program/src/bedtools/bedtools2/profile
rpm=$(zcat ${DPATH}/${File}_1.fastq.gz|wc -l|cut -f 1 -d ' '|awk '{printf "%f\n",$1/4000000.0}')
rpm=$(awk '{printf "%f\n",1/("'$rpm'")}' /sibcb1/wuweilab1/wuwei2022/Wuwei/Pipeline/tst.txt)

bedtools genomecov -bg -strand "+" -scale $rpm -du -ibam $DPATH/$File.$SZ.rmDup.sorted.bam|bedtools intersect -a - -b $BLACKLIST -v |bedtools sort -i stdin> $DPATH/$File.$SZ.pos.bedgraph
bedtools genomecov -bg -strand "-" -scale $rpm -du -ibam $DPATH/$File.$SZ.rmDup.sorted.bam|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools intersect -a - -b $BLACKLIST -v |bedtools sort -i stdin> $DPATH/$File.$SZ.neg.bedgraph
bedGraphToBigWig $DPATH/$File.$SZ.pos.bedgraph $CHROMSIZE $DPATH/$File.$SZ.pos.bw
bedGraphToBigWig $DPATH/$File.$SZ.neg.bedgraph $CHROMSIZE $DPATH/$File.$SZ.neg.bw
rm $DPATH/$File.$SZ.pos.bedgraph
rm $DPATH/$File.$SZ.neg.bedgraph
