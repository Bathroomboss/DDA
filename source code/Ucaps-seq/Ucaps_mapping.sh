#!/bin/bash
# Ucaps_mapping



cutadapt -j 24 --discard-trimmed -m 55 -g ^GACTGGTTCCAATTGAAAGTGCTCTTCCGATCT -o ${sample}_cutadapt_1.fastq.gz -p ${sample}_cutadapt_2.fastq.gz ${sample}_1.fastq.gz ${sample}_2.fastq.gz
bowtie -m 1 -p 24 $BOWTIEINDEX -1 ${sample}_cutadapt_1.fastq.gz -2 ${sample}_cutadapt_2.fastq.gz -S ${sample}.${SZ}.sam
samtools view -f 67 -q 60 -bSh -@ 24 ${sample}.${SZ}.sam > ${sample}.${SZ}.bam 
rm ${sample}.${SZ}.sam
samtools sort -@ 24 -o ${sample}.${SZ}.sorted.bam ${sample}.${SZ}.bam
rm ${sample}.${SZ}.bam
java -Xms20g -Xmx20g -XX:ParallelGCThreads=4 -jar /sibcb1/wuweilab1/liangyu/Software/picard/picard.jar MarkDuplicates I=${sample}.${SZ}.sorted.bam O=${sample}.${SZ}.rmDup.sorted.bam M=${sample}.rmDup_metrics.txt ASO=coordinate REMOVE_DUPLICATES=true 
rm ${sample}.${SZ}.sorted.bam
bedtools bamtobed -i ${sample}.${SZ}.rmDup.sorted.bam |bedtools intersect -v -a - -b ${BLACKLIST} > ${sample}.${SZ}.sorted.bed 
awk 'BEGIN{FS=OFS="\t"}{if ($6=="+" && $2>=1) {print $1,$2-1,$2,$4,$5,$6}else if ($6=="-"){print $1,$3,$3+1,$4,$5,$6}}' ${sample}.$SZ.sorted.bed > ${sample}.$SZ.breaksite.sorted.bed 
rpm=$(wc -l ${sample}.$SZ.breaksite.sorted.bed |cut -f 1 -d ' '|awk '{printf "%f\n",$1/1000000.0}')
rpm=$(awk '{printf "%f\n",1/("'$rpm'")}' /sibcb1/wuweilab1/wuwei2022/Wuwei/Pipeline/tst.txt)
bedtools genomecov -bg -strand "+" -scale $rpm -i ${sample}.$SZ.breaksite.sorted.bed -g $CHROMSIZE|bedtools sort -i stdin> ${sample}.${SZ}.pos.bedgraph
bedtools genomecov -bg -strand "-" -scale $rpm -i ${sample}.$SZ.breaksite.sorted.bed -g $CHROMSIZE|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin > ${sample}.${SZ}.neg.bedgraph
bedGraphToBigWig ${sample}.${SZ}.pos.bedgraph $CHROMSIZE ${sample}.${SZ}.pos.bw
bedGraphToBigWig ${sample}.${SZ}.neg.bedgraph $CHROMSIZE ${sample}.${SZ}.neg.bw
rm ${sample}.${SZ}.pos.bedgraph ${sample}.${SZ}.neg.bedgraph

awk 'BEGIN{FS=OFS="\t"}{if ($2 >=5) print $1,$2-5,$3+5,$4,$5,$6}' ${sample}.$SZ.breaksite.sorted.bed > ${sample}.$SZ.breaksite_5bp.sorted.bed
bedtools getfasta -fi $FA -bed ${sample}.$SZ.breaksite_5bp.sorted.bed -fo ${sample}.$SZ.breaksite_5bp.fa -s 
awk '{if (NR%2==0) print $0}' ${sample}.$SZ.breaksite_5bp.fa > ${sample}.$SZ.breaksite_5bp.seq.txt
Rscript /sibcb2/wuweilab2/liangyu/Pipeline/Ucaps.r ${sample}.$SZ


date  