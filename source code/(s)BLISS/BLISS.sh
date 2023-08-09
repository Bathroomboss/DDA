#!/bin/bash
#(s)BLISS
# Note: 
# The information of barcode is important but cannot be determined from the data of some samples.
# And for these samples, we just use bowtie2 --local to perform analysis.
# Otherwise, we extract umi and remove duplicates via umi_tools

umi_tools extract --stdin=$File.fastq.gz --stdout=$File.processed.fastq.gz \
--extract-method=regex --bc-pattern="(?P<umi_1>.{8})(?P<discard_1>$Barcode{s<=1})" \
--log=$File_processed.log

bowtie --best -p 24 -n 3 -l 50 -k 1 $BowtieIndex -q $File.processed.fastq.gz -S $File.sam --chunkmbs 200 --no-unal
samtools view -bSh $File.sam  > $File.bam
rm ${File}.sam
samtools sort -@ 24 -o $File.sorted.bam $File.bam
rm ${File}.bam
samtools index $File.sorted.bam
umi_tools dedup -I $File.sorted.bam -S $File.dedup.sorted.bam -L $File.dedup.log --mapping-quality=30
bedtools bamtobed -i ${File}.dedup.sorted.bam|bedtools intersect -v -a - -b ${Blacklist} > ${File}.sorted.bed

rpm=$(wc -l ${File}.sorted.bed|cut -f 1 -d ' '|awk '{printf "%f\n",$1/1000000.0}')
rpm=$(awk '{printf "%f\n",1/("'$rpm'")}' /sibcb1/wuweilab1/wuwei2022/Wuwei/Pipeline/tst.txt)

bedtools genomecov -bg -strand "+" -scale $rpm -i ${File}.sorted.bed -g ${chromSize}|bedtools sort -i stdin> ${File}.pos.bedgraph
bedtools genomecov -bg -strand "-" -scale $rpm -i ${File}.sorted.bed -g ${chromSize}|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin > ${File}.neg.bedgraph
bedtools genomecov -bg -5 -strand "+" -scale $rpm -i ${File}.sorted.bed -g ${chromSize}|bedtools sort -i stdin> ${File}.pos5.bedgraph
bedtools genomecov -bg -5 -strand "-" -scale $rpm -i ${File}.sorted.bed -g ${chromSize}|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin > ${File}.neg5.bedgraph

bedGraphToBigWig ${File}.pos.bedgraph ${chromSize} ${File}.pos.bw
rm ${File}.pos.bedgraph
bedGraphToBigWig ${File}.neg.bedgraph ${chromSize} ${File}.neg.bw
rm ${File}.neg.bedgraph
bedGraphToBigWig ${File}.pos5.bedgraph ${chromSize} ${File}.pos5.bw
rm ${File}.pos5.bedgraph
bedGraphToBigWig ${File}.neg5.bedgraph ${chromSize} ${File}.neg5.bw
rm ${File}.neg5.bedgraph

date
