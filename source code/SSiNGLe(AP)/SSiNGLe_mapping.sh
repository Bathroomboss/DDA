#!/bin/bash
# SSiNGLe

cutadapt -j 0 -g ^AGTTGCGGATGGGGGGGGGG -G ^TTTTTTTTTTTT \
-o ${File}_cuthead_1.fastq.gz  -p ${File}_cuthead_2.fastq.gz ${File}_1.fastq.gz ${File}_2.fastq.gz --discard-untrimmed
cutadapt -j 0 -a AGATCGGAAGAG -A AGATCGGAAGAG -q 20 \
-o ${File}_trim_1.fastq.gz -p ${File}_trim_2.fastq.gz ${File}_cuthead_1.fastq.gz ${File}_cuthead_2.fastq.gz

rm  ${File}_cuthead_1.fastq.gz ${File}_cuthead_2.fastq.gz

bwa mem -t 24 -M -R '@RG\tID:foo\tSM:bar\tLB:library1' $bwa_Index ${File}_trim_1.fastq.gz ${File}_trim_2.fastq.gz > ${File}.$SZ.sam
grep '@' ${File}.$SZ.sam > ${File}.$SZ.head.txt
awk '{if ($2 == 147) print $0 ;else if ($2 == 163) print $0}' ${File}.$SZ.sam |awk '{if ($5 > 20) print $0}' > ${File}.filter.txt
cat ${File}.$SZ.head.txt ${File}.filter.txt > ${File}.filter.sam
samtools view -@ 24 -bS ${File}.filter.sam > ${File}.$SZ.bam
rm ${File}.$SZ.head.txt ${File}.filter.txt
rm ${File}.filter.sam
samtools view -@ 24 -bSh ${File}.$SZ.sam > ${File}.$SZ.raw.bam
rm ${File}.$SZ.sam

# extract only unique and proper mapped paired reads and extract reads2
samtools sort -@ 24 -o ${File}.$SZ.sorted.bam ${File}.$SZ.bam
rm ${File}.$SZ.bam
bedtools bamtobed -tag NM -i ${File}.$SZ.sorted.bam |bedtools intersect -v -a - -b ${Blacklist} > ${File}.sorted.bed
# we use the -tag option to select the BAM edit distance (the NM tag) as the score column in the resulting BED records.

awk 'BEGIN{FS=OFS="\t"}{if ($6 == "+") print $1,$2,$2+1,$4,$5,"-";else print $1,$3-1,$3,$4,$5,"+"}' ${File}.sorted.bed | awk '$6 != "-" || $2 > 20 {print}'> ${File}.breaksite.bed
awk 'BEGIN{FS=OFS="\t"}{if ($6 == "+") print $1,$3,$3+20,$4,$5,$6;else print $1,$2-20,$2,$4,$5,$6}' ${File}.breaksite.bed > ${File}.20bp_downstream.bed

bedtools getfasta -fi $fa -bed ${File}.20bp_downstream.bed -fo ${File}.20bp_downstream.fa
python /sibcb2/wuweilab2/liangyu/SSiNGLe/SSiNGLe_retain.py ${File}.20bp_downstream.fa ${File}.20bp_downstream_retain.txt
awk 'NR%2==1' ${File}.20bp_downstream_retain.txt > ${File}.20bp_downstream_retain.bed
rm ${File}.20bp_downstream_retain.txt
paste ${File}.20bp_downstream.bed ${File}.breaksite.bed |bedtools intersect -a - -b ${File}.20bp_downstream_retain.bed -u |awk 'BEGIN{FS=OFS="\t"}{print $7,$8,$9,$10,$11,$12}' > ${File}.$SZ.filtered.breaksite.bed
rm ${File}.20bp_downstream.bed ${File}.breaksite.bed ${File}.20bp_downstream.fa ${File}.20bp_downstream_retain.bed
sort -k1,1 -k2,2n ${File}.$SZ.filtered.breaksite.bed > ${File}.$SZ.filtered.breaksite.sorted.bed
rm ${File}.$SZ.filtered.breaksite.bed
rpm=$(wc -l $File.$SZ.filtered.breaksite.sorted.bed|cut -f 1 -d ' '|awk '{printf "%f\n",$1/1000000.0}')
rpm=$(awk '{printf "%f\n",1/("'$rpm'")}' /sibcb1/wuweilab1/wuwei2022/Wuwei/Pipeline/tst.txt)
bedtools genomecov -bg -strand "+" -scale $rpm -i ${File}.$SZ.filtered.breaksite.sorted.bed -g ${chromSize}|bedtools sort -i stdin> ${File}.pos.bedgraph
bedtools genomecov -bg -strand "-" -scale $rpm -i ${File}.$SZ.filtered.breaksite.sorted.bed -g ${chromSize}|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin > ${File}.neg.bedgraph
bedGraphToBigWig ${File}.pos.bedgraph ${chromSize} ${File}.pos.bw
rm ${File}.pos.bedgraph
bedGraphToBigWig ${File}.neg.bedgraph ${chromSize} ${File}.neg.bw
rm ${File}.neg.bedgraph

date




