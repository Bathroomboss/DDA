#!/bin/bash
# GLOE-seq + Endoseq

date
if [ ${END}=="SE" ] ; then
        java -jar trimmomatic-0.39.jar SE -threads 24 \
        ${File}.fastq.gz ${File}.trim.fastq.gz \
        ILLUMINACLIP:adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        bowtie2 -p 24 -x $Bowtie2Index -U ${File}.trim.fastq.gz -S $File.sam
        samtools view -bSh -@ 24 -q 30 $File.sam > $File.bam
else
        java -jar trimmomatic-0.39.jar PE -threads 24 \
        ${File}_1.fastq.gz ${File}_2.fastq.gz \
        ${File}_trim_1.fastq.gz ${File}_trim_unpaired_1.fastq.gz ${File}_trim_2.fastq.gz ${File}_trim_unpaired_2.fastq.gz \
        ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        bowtie2 -p 24 -x $Bowtie2Index -1 ${File}_trim_1.fastq.gz -2 ${File}_trim_2.fastq.gz -S ${File}.sam
        samtools view -bSh -@ 24 -q 30 -f 64 ${File}.sam > $File.bam
fi
samtools sort -@ 24 -o $File.$SZ.sorted.bam $File.bam
bamToBed -i $File.sorted.bam > $File.$SZ.raw.bed
awk 'BEGIN{FS=OFS="\t"}{if ($6=="+") print $1,$2,$3,$4,$5,"-" ;else print $1,$2,$3,$4,$5,"+"}'  $File.$SZ.raw.bed > $File.$SZ.sorted.bed
rm $File.sam $File.bam

rpm=$(wc -l $File.$SZ.sorted.bed|cut -f 1 -d ' '|awk '{printf "%f\n",$1/1000000.0}')
rpm=$(awk '{printf "%f\n",1/("'$rpm'")}' /sibcb1/wuweilab1/wuwei2022/Wuwei/Pipeline/tst.txt)

bedtools genomecov -bg -strand "+" -scale $rpm -i ${File}.$SZ.sorted.bed -g ${chromSize}|bedtools sort -i stdin> ${File}.pos.bedgraph
bedtools genomecov -bg -strand "-" -scale $rpm -i ${File}.$SZ.sorted.bed -g ${chromSize}|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin > ${File}.neg.bedgraph
bedtools genomecov -bg -3 -strand "+" -scale $rpm -i ${File}.$SZ.sorted.bed -g ${chromSize}|bedtools sort -i stdin> ${File}.pos3.bedgraph
bedtools genomecov -bg -3 -strand "-" -scale $rpm -i ${File}.$SZ.sorted.bed -g ${chromSize}|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin > ${File}.neg3.bedgraph

bedGraphToBigWig ${File}.pos.bedgraph ${chromSize} ${File}.pos.bw
rm ${File}.pos.bedgraph
bedGraphToBigWig ${File}.neg.bedgraph ${chromSize} ${File}.neg.bw
rm ${File}.neg.bedgraph
bedGraphToBigWig ${File}.pos3.bedgraph ${chromSize} ${File}.pos3.bw
rm ${File}.pos3.bedgraph
bedGraphToBigWig ${File}.neg3.bedgraph ${chromSize} ${File}.neg3.bw
rm ${File}.neg3.bedgraph