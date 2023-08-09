#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -j y
#$ -m n
#$ -M wuw@sibcb.ac.cn
#$ -notify
#$ -N HydEn_mapping

export PATH=$PATH:/sibcb1/wuweilab1/liangyu/miniconda3/bin
source activate


cutadapt -j 24 -m 15 -q 10 --match-read-wildcards -a AGATCGGAAGAG -A AGATCGGAAGAG -o ${sample}_cutadapt_1.fastq.gz -p ${sample}_cutadapt_2.fastq.gz ${sample}_1.fastq.gz ${sample}_2.fastq.gz
bowtie -p 24 -v 1 -m 2 /sibcb1/wuweilab1/wuwei-lab/Data/DamageView/RawData/HydEn-seq/GSE62181/Oligo_BTIndex/Oligo -q ${sample}_cutadapt_1.fastq.gz --un ${sample}_cutadapt_unalinOligo_1.fastq -S ${sample}_cutadapt_1_alinOligo.sam
rm ${sample}_cutadapt_1_alinOligo.sam
gunzip ${sample}_cutadapt_2.fastq.gz 
fastq_pair ${sample}_cutadapt_unalinOligo_1.fastq ${sample}_cutadapt_2.fastq
bowtie -v 1 -m 2 -p 24 -X 10000 --best ${L03_BTIndex} -1 ${sample}_cutadapt_unalinOligo_1.fastq.paired.fq -2 ${sample}_cutadapt_2.fastq.paired.fq -S ${sample}.${SZ}_PE.sam

rm ${sample}_cutadapt_unalinOligo_1.fastq.single.fq ${sample}_cutadapt_2.fastq.single.fq 
gzip -1  ${sample}_cutadapt_unalinOligo_1.fastq.paired.fq ${sample}_cutadapt_2.fastq.paired.fq 
gzip -1 ${sample}_cutadapt_2.fastq ${sample}_cutadapt_unalinOligo_1.fastq 

samtools view -@ 24 -f 4 -bSh ${sample}.${SZ}_PE.sam > ${sample}.${SZ}_PEunmap.bam
samtools fastq -@ 24 -f 64 ${sample}.${SZ}_PEunmap.bam > ${sample}.${SZ}_PEunmap_1.fastq
bowtie -v 1 -m 2 -p 24 ${L03_BTIndex} -q ${sample}.${SZ}_PEunmap_1.fastq -S ${sample}.${SZ}_PEunmap_SEmap.sam 
gzip -1 ${sample}.${SZ}_PEunmap_1.fastq

samtools view -f 67 -@ 24 -bSh ${sample}.${SZ}_PE.sam > ${sample}.${SZ}_PE_R1.bam
samtools view -@ 24 -bSh ${sample}.${SZ}_PEunmap_SEmap.sam > ${sample}.${SZ}_PEunmap_SEmap.bam 
rm ${sample}.${SZ}_PE.sam  ${sample}.${SZ}_PEunmap_SEmap.sam  

samtools sort -@ 24 -o ${sample}.${SZ}_PE_R1.sorted.bam ${sample}.${SZ}_PE_R1.bam
samtools sort -@ 24 -o ${sample}.${SZ}_PEunmap_SEmap.sorted.bam ${sample}.${SZ}_PEunmap_SEmap.bam

rm ${sample}.${SZ}_PE_R1.bam ${sample}.${SZ}_PEunmap_SEmap.bam

bedtools bamtobed -i ${sample}.${SZ}_PE_R1.sorted.bam|bedtools intersect -v -a - -b ${BLACKLIST} > ${sample}.${SZ}_PE_R1.sorted.bed
bedtools bamtobed -i ${sample}.${SZ}_PEunmap_SEmap.sorted.bam|bedtools intersect -v -a - -b ${BLACKLIST} > ${sample}.${SZ}_PEunmap_SEmap.sorted.bed

cat ${sample}.${SZ}_PE_R1.sorted.bed ${sample}.${SZ}_PEunmap_SEmap.sorted.bed |bedtools sort -i - |awk 'BEGIN{FS=OFS="\t"}{if ($6=="+" && $2>=1) {print $1,$2-1,$3,$4,$5,$6}else if ($6=="-"){print $1,$2,$3+1,$4,$5,$6}}' > ${sample}.${SZ}.sorted.bed

rpm=$(wc -l ${sample}.${SZ}.sorted.bed|cut -f 1 -d ' '|awk '{printf "%f\n",$1/1000000.0}')
rpm=$(awk '{printf "%f\n",1/("'$rpm'")}' /sibcb1/wuweilab1/wuwei2022/Wuwei/Pipeline/tst.txt)
bedtools genomecov -bg -strand "+" -scale $rpm -i ${sample}.${SZ}.sorted.bed -g $CHROMSIZE|bedtools sort -i stdin> ${sample}.${SZ}.pos.bedgraph
bedtools genomecov -bg -strand "-" -scale $rpm -i ${sample}.${SZ}.sorted.bed -g $CHROMSIZE|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin > ${sample}.${SZ}.neg.bedgraph

bedtools genomecov -5 -bg -strand "+" -scale $rpm -i ${sample}.${SZ}.sorted.bed -g $CHROMSIZE|bedtools sort -i stdin> ${sample}.${SZ}.pos5.bedgraph
bedtools genomecov -5 -bg -strand "-" -scale $rpm -i ${sample}.${SZ}.sorted.bed -g $CHROMSIZE|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin > ${sample}.${SZ}.neg5.bedgraph

bedGraphToBigWig ${sample}.${SZ}.pos.bedgraph $CHROMSIZE ${sample}.${SZ}.pos.bw
bedGraphToBigWig ${sample}.${SZ}.pos5.bedgraph $CHROMSIZE ${sample}.${SZ}.pos5.bw
bedGraphToBigWig ${sample}.${SZ}.neg.bedgraph $CHROMSIZE ${sample}.${SZ}.neg.bw
bedGraphToBigWig ${sample}.${SZ}.neg5.bedgraph $CHROMSIZE ${sample}.${SZ}.neg5.bw

rm ${sample}.${SZ}.pos.bedgraph ${sample}.${SZ}.pos5.bedgraph ${sample}.${SZ}.neg.bedgraph ${sample}.${SZ}.neg5.bedgraph

bedtools genomecov -bg -scale $rpm -i ${sample}.${SZ}.sorted.bed -g $CHROMSIZE|bedtools sort -i stdin> ${sample}.${SZ}.bedgraph
bedGraphToBigWig ${sample}.${SZ}.bedgraph $CHROMSIZE ${sample}.${SZ}.bw
rm ${sample}.${SZ}.bedgraph
date