#!/bin/bash
# TrAEL-seq

umi_tools extract --stdin=$File.fastq.gz --stdout=$File.markUMI.fastq.gz --extract-method=regex --bc-pattern='(?P<umi_1>.{8})'

gunzip $File.markUMI.fastq.gz
python TrAEL_process.py $File.markUMI.fastq $File.trimT.fastq
rm $File.markUMI.fastq
gzip $File.trimT.fastq
bowtie2 -p 24 -x ${Bowtie2Index} -U $File.trimT.fastq.gz -S $File.$SZ.sam --local
samtools view -@ 24 -bSh $File.$SZ.sam > $File.$SZ.bam
rm ${File}.$SZ.sam
samtools sort -@ 24 -o $File.$SZ.sorted.bam $File.$SZ.bam
rm ${File}.$SZ.bam
samtools index $File.$SZ.sorted.bam
umi_tools dedup -I $File.$SZ.sorted.bam -S $File.$SZ.dedup.sorted.bam
bedtools bamtobed -i ${File}.$SZ.dedup.sorted.bam > ${File}.$SZ.dedup.sorted.temp.bed
awk 'BEGIN{FS=OFS="\t"}{if ($6=="+") print $1,$2,$3,$4,$5,"-" ;else print $1,$2,$3,$4,$5,"+"}' ${File}.$SZ.dedup.sorted.temp.bed > ${File}.$SZ.dedup.sorted.bed
intersectBed -a ${File}.$SZ.dedup.sorted.bed -b ${Blacklist} -v > ${File}.$SZ.dedup.sorted.rmblacklist.bed
rm ${File}.$SZ.dedup.sorted.bed ${File}.$SZ.dedup.sorted.temp.bed

rpm=$(wc -l ${File}.$SZ.dedup.sorted.rmblacklist.bed|cut -f 1 -d ' '|awk '{printf "%f\n",$1/1000000.0}')
rpm=$(awk '{printf "%f\n",1/("'$rpm'")}' /sibcb1/wuweilab1/wuwei2022/Wuwei/Pipeline/tst.txt)

bedtools genomecov -bg -strand "+" -scale $rpm -i ${File}.$SZ.dedup.sorted.rmblacklist.bed -g ${chromSize}|bedtools sort -i stdin> ${File}.pos.bedgraph
bedtools genomecov -bg -strand "-" -scale $rpm -i ${File}.$SZ.dedup.sorted.rmblacklist.bed -g ${chromSize}|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin > ${File}.neg.bedgraph
bedtools genomecov -bg -3 -strand "+" -scale $rpm -i ${File}.$SZ.dedup.sorted.rmblacklist.bed -g ${chromSize}|bedtools sort -i stdin> ${File}.pos3.bedgraph
bedtools genomecov -bg -3 -strand "-" -scale $rpm -i ${File}.$SZ.dedup.sorted.rmblacklist.bed -g ${chromSize}|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin > ${File}.neg3.bedgraph

bedGraphToBigWig ${File}.pos.bedgraph ${chromSize} ${File}.pos.bw
rm ${File}.pos.bedgraph
bedGraphToBigWig ${File}.neg.bedgraph ${chromSize} ${File}.neg.bw
rm ${File}.neg.bedgraph
bedGraphToBigWig ${File}.pos3.bedgraph ${chromSize} ${File}.pos3.bw
rm ${File}.pos3.bedgraph
bedGraphToBigWig ${File}.neg3.bedgraph ${chromSize} ${File}.neg3.bw
rm ${File}.neg3.bedgraph

date




