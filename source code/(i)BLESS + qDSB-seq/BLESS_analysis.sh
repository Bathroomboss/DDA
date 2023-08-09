#!/bin/bash
# BLESS
# Note: 
# The information of barcode is important but cannot be determined from the data of some samples.
# And for these samples, we just use bowtie2 --local to perform analysis.
# Otherwise, we cut the barcode part viy cutadapt firstly and use bowtie with parameters of "--best -p 24  -m 1 -v 1"    
# The following code is an example and it should be noted that we did not remove duplicates for single-ended sequenced files
date
cutadapt -e 0.1 -O 3 -m 23 --discard-untrimmed --quality-cutoff 20 -g TCGAGACGACG -g TCGAGGTAGTA -G TCGAGACGACG -G TCGAGGTAGTA  \
-o ${Filename}_cutadapt_1.fastq.gz -p ${Filename}_cutadapt_2.fastq.gz ${Filename}_1.fastq.gz ${Filename}_2.fastq.gz

bowtie --best -p 24  -m 1 -v 1  ${BowtieIndex} -1 ${Filename}_cutadapt_1.fastq.gz -2 ${Filename}_cutadapt_2.fastq.gz -S ${File}.sam --chunkmbs 200 --no-unal
samtools view -bhS ${File}.sam > ${File}.bam
samtools sort -@ 24 -o  ${File}.sorted.bam ${File}.bam

java -jar -Xms20g -Xmx20g -XX:ParallelGCThreads=4  /sibcb/program/install/picard/picard.jar  \
MarkDuplicates I=${File}.sorted.bam  \
O=${File}.sorted.rmdup.bam M=${File}.rmdup_matrix.txt REMOVE_DUPLICATES=true

samtools view -h -f 64 ${File}.sorted.rmdup.bam|samtools view -bS - > ${File}.sorted.rmdup.R1.bam
bedtools bamtobed -i ${File}.sorted.rmdup.R1.bam|bedtools intersect -v -a - -b ${Blacklist} > ${File}.sorted.bed

rpm=$(wc -l ${File}.sorted.bed|cut -f 1 -d ' '|awk '{printf "%f\n",$1/1000000.0}')
rpm=$(awk '{printf "%f\n",1/("'$rpm'")}' /sibcb1/wuweilab1/wuwei2022/Wuwei/Pipeline/tst.txt)

bedtools genomecov -bg -strand "+" -scale $rpm -i ${File}.sorted.bed -g ${chromSize}|bedtools sort -i stdin> ${File}.pos.bedgraph
bedtools genomecov -bg -strand "-" -scale $rpm -i ${File}.sorted.bed -g ${chromSize}|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin > ${File}.neg.bedgra
ph
bedtools genomecov -bg -5 -strand "+" -scale $rpm -i ${File}.sorted.bed -g ${chromSize}|bedtools sort -i stdin> ${File}.pos5.bedgraph
bedtools genomecov -bg -5 -strand "-" -scale $rpm -i ${File}.sorted.bed -g ${chromSize}|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin > ${File}.neg5.be
dgraph

bedGraphToBigWig ${File}.pos.bedgraph ${chromSize} ${File}.pos.bw
bedGraphToBigWig ${File}.neg.bedgraph ${chromSize} ${File}.neg.bw
bedGraphToBigWig ${File}.pos5.bedgraph ${chromSize} ${File}.pos5.bw
bedGraphToBigWig ${File}.neg5.bedgraph ${chromSize} ${File}.neg5.bw

rm ${File}.bam
rm ${File}.sam
rm ${File}.pos.bedgraph
rm ${File}.neg.bedgraph
rm ${File}.pos5.bedgraph
rm ${File}.neg5.bedgraph



