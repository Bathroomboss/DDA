#!/bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -j y
#$ -m n
#$ -M wuw@sibcb.ac.cn
#$ -notify
#$ -N Damage-seq_mapping

date
BOWTIEINDEX="/sibcb1/wuweilab1/wuwei-lab/Reference/Human/hg19/BowtieIndexplusrDNA/genome"
BLACKLIST="/sibcb1/wuweilab1/wuwei2022/Wuwei/Pipeline/hg19-blacklist.v2.bed"
CHROMSIZE="/sibcb1/wuweilab1/wuwei-lab/Reference/Human/hg19/hg19plusrDNA.chrom.sizes"
SZ=hg19
FA=/sibcb1/wuweilab1/wuwei-lab/Reference/Human/hg19/hg19plusrDNA.fa

source /sibcb/program/install/python-3.9/profile
python3 -m pip install --user cutadapt
cutadapt -j 24 --discard-trimmed -m 15 -g ^GACTGGTTCCAATTGAAAGTGCTCTTCCGATCT  -o ${sample}_cutadapt_1.fastq.gz -p ${sample}_cutadapt_2.fastq.gz ${sample}_1.fastq.gz ${sample}_2.fastq.gz

bowtie -p 24 -X 1000 -m 4 --seed 123 --nomaqround ${BOWTIEINDEX} -1 ${sample}_cutadapt_1.fastq.gz -2 ${sample}_cutadapt_2.fastq.gz -S ${sample}.$SZ.sam --chunkmbs 200 
samtools view -@ 24 -F 1024 -bSh ${sample}.$SZ.sam > ${sample}.$SZ.bam
rm ${sample}.$SZ.sam 
samtools sort -@ 24 -o ${sample}.$SZ.sorted.bam ${sample}.$SZ.bam
rm ${sample}.$SZ.bam
bedtools bamtobed -i ${sample}.$SZ.sorted.bam  |bedtools intersect -v -a - -b $BLACKLIST |grep -v "\/2" > ${sample}.$SZ.sorted.bed  
awk 'BEGIN{FS=OFS="\t"}{if ($6=="+" && $2>=2) {print $1,$2-2,$2,$4,$5,$6}else if ($6=="-"){print $1,$3,$3+2,$4,$5,$6}}' ${sample}.$SZ.sorted.bed > ${sample}.$SZ.breaksite.sorted.bed  
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
source /sibcb/program/install/r-4.1/profile
Rscript /sibcb2/wuweilab2/liangyu/Pipeline/Damageseq.r ${sample}.$SZ 

date 


