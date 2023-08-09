#!/bin/bash
# S1-seq + ExoT-seq


echo $File

if [ ${END} == "SE" ];then
        bowtie2 -p 24 -x $Bowtie2Index -N 1 -X 1000 -U ${DPATH}/${File}.fastq.gz -S ${DPATH}/${File}.$SZ.sam --no-unal
else
        bowtie2 -p 24 -x $Bowtie2Index -N 1 -X 1000 -1 ${DPATH}/${File}_1.fastq.gz -2 ${DPATH}/${File}_2.fastq.gz -S ${DPATH}/${File}.$SZ.sam --no-unal
fi

samtools view -@ 24 -F 2820 -bSh ${DPATH}/${File}.${SZ}.sam  > ${DPATH}/${File}.${SZ}.bam
rm ${DPATH}/${File}.${SZ}.sam
samtools sort -@ 24 -o ${DPATH}/${File}.${SZ}.sorted.bam ${DPATH}/${File}.${SZ}.bam
rm ${File}.${SZ}.bam
if [ ${END} == "SE" ] ;then
        bedtools bamtobed -i ${DPATH}/${File}.${SZ}.sorted.bam|bedtools intersect -v -a - -b $Blacklist  > ${DPATH}/${File}.${SZ}.sorted.bed
else
        java -Xms20g -Xmx20g -XX:ParallelGCThreads=8 -jar /sibcb/program/install/picard/picard.jar MarkDuplicates I=${DPATH}/${File}.${SZ}.sorted.bam O=${DPATH}/${File}.${SZ}.rmDup.sorted.bam M=${DPATH}/${File}.rmDup_metrics.txt ASO=coordinate REMOVE_DUPLICATES=true
        samtools view -f 65 -@ 24 -q 10 -F 1024 -bSh ${File}.${SZ}.rmDup.sorted.bam >  ${File}.${SZ}.rmDup.sorted.R1.bam
        bedtools bamtobed -i ${File}.${SZ}.rmDup.sorted.R1.bam|bedtools intersect -v -a - -b $Blacklist  > ${File}.${SZ}.sorted.bed
fi
rpm=$(wc -l $DPATH/$File.$SZ.sorted.bed |cut -f 1 -d ' '|awk '{printf "%f\n",$1/1000000.0}')
rpm=$(awk '{printf "%f\n",1/("'$rpm'")}' /sibcb1/wuweilab1/wuwei2022/Wuwei/Pipeline/tst.txt)
bedtools genomecov -bg -strand "+" -scale $rpm -i $DPATH/$File.$SZ.sorted.bed -g $chromSize|bedtools sort -i stdin> ${DPATH}/${File}.pos.bedgraph
bedtools genomecov -bg -strand "-" -scale $rpm -i $DPATH/$File.$SZ.sorted.bed -g $chromSize|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin > ${DPATH}/${File}.neg.bedgraph
bedtools genomecov -bg -5 -strand "+" -scale $rpm -i $DPATH/$File.$SZ.sorted.bed -g $chromSize|bedtools sort -i stdin> ${DPATH}/${File}.pos5.bedgraph
bedtools genomecov -bg -5 -strand "-" -scale $rpm -i $DPATH/$File.$SZ.sorted.bed -g $chromSize|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin > ${DPATH}/${File}.neg5.bedgraph
bedGraphToBigWig ${DPATH}/${File}.pos.bedgraph $chromSize ${DPATH}/${File}.$SZ.pos.bw
bedGraphToBigWig ${DPATH}/${File}.neg.bedgraph $chromSize ${DPATH}/${File}.$SZ.neg.bw
bedGraphToBigWig ${DPATH}/${File}.pos5.bedgraph $chromSize ${DPATH}/${File}.$SZ.pos5.bw
bedGraphToBigWig ${DPATH}/${File}.neg5.bedgraph $chromSize ${DPATH}/${File}.$SZ.neg5.bw
rm ${DPATH}/${File}.pos.bedgraph
rm ${DPATH}/${File}.neg.bedgraph
rm ${DPATH}/${File}.pos5.bedgraph
rm ${DPATH}/${File}.neg5.bedgraph
