#!/bin/bash
# Spo11-oligo_mapping

# The following code for cutting oligos introduced by the library prepartion applies to most data. But we recommend you use fastqc to look at the sequence characteristics at the first.

# Remove the first 5 NTs from R1, as well as the portion after GGGAGAT.
cutadapt -j 24 -u 5 -a GGGAGAT -o ${sample}_Trim1_1.fastq.gz -p ${sample}_Trim1_2.fastq.gz ${sample}_1.fastq.gz ${sample}_2.fastq.gz
# Remove the leading CCC from R1 and the portion after GGGNNNNNAGAT from R2(including GGGNNNNNAGAT).
cutadapt -j 24 -g ^CCC -A GGGNNNNNAGAT --discard-untrimmed -o ${sample}_Trim2_1.fastq.gz -p ${sample}_Trim2_2.fastq.gz ${sample}_Trim1_1.fastq.gz ${sample}_Trim1_2.fastq.gz
# Remove the leading CCC from R2.
cutadapt -j 24 -G ^CCC -o ${sample}_Trim3_1.fastq.gz -p ${sample}_Trim3_2.fastq.gz ${sample}_Trim2_1.fastq.gz ${sample}_Trim2_2.fastq.gz
# Remove the leading G (if present) from R1 and R2.
cutadapt -e 0 -j 24 -g ^C -G ^C -m 15 -o ${sample}_Trim4_1.fastq.gz -p ${sample}_Trim4_2.fastq.gz ${sample}_Trim3_1.fastq.gz ${sample}_Trim3_2.fastq.gz 
# Remove the trailing G (if present) from R1 and R2.
cutadapt -e 0 -j 24 -a G$ -A G$ -m 15 -o ${sample}_Trim5_1.fastq.gz -p ${sample}_Trim5_2.fastq.gz ${sample}_Trim4_1.fastq.gz ${sample}_Trim4_2.fastq.gz

rm ${sample}_Trim1_1.fastq.gz ${sample}_Trim1_2.fastq.gz 
rm ${sample}_Trim2_1.fastq.gz ${sample}_Trim2_2.fastq.gz
rm ${sample}_Trim3_1.fastq.gz ${sample}_Trim3_2.fastq.gz 
rm ${sample}_Trim4_1.fastq.gz ${sample}_Trim4_2.fastq.gz

#################### For SE cutadapt ################### #################### For SE cutadapt ###################
cutadapt -j 24 -u 5 -a GGGAGAT -o ${sample}_Trim1.fastq.gz ${sample}.fastq.gz
cutadapt -j 24 -g ^CCC --discard-untrimmed -o ${sample}_Trim2.fastq.gz ${sample}_Trim1.fastq.gz
cutadapt -e 0 -j 24 -g ^C -o ${sample}_Trim3.fastq.gz ${sample}_Trim2.fastq.gz 
cutadapt -e 0 -j 24 -a G$ -m 15 -o ${sample}_Trim4.fastq.gz ${sample}_Trim3.fastq.gz
rm ${sample}_Trim1.fastq.gz ${sample}_Trim2.fastq.gz ${sample}_Trim3.fastq.gz
#################### For SE cutadapt ################### #################### For SE cutadapt ###################

# Note: we perform both "--end-to-end" and "--local" model in bowtie2 and choose the results with higher mapping rate. 
bowtie2 -p 24 -x $Bowtie2Index -X 1000 --no-discordant --very-sensitive --mp 5,1 --np 0 -1 ${sample}_Trim5_1.fastq.gz -2 ${sample}_Trim5_2.fastq.gz -S ${sample}.mm10.sam --no-unal
bowtie2 -p 24 -x $Bowtie2Index -X 1000 --no-discordant --very-sensitive --mp 5,1 --np 0 -1 ${sample}_Trim5_1.fastq.gz -2 ${sample}_Trim5_2.fastq.gz -S ${sample}.mm10_local.sam --local --no-unal

samtools view -bSh -@ 24 ${sample}.${SZ}.sam > ${sample}.${SZ}.bam
samtools sort -@ 24 -o ${sample}.${SZ}.sorted.bam ${sample}.${SZ}.bam
rm ${sample}.${SZ}.bam
bedtools bamtobed -i ${sample}.${SZ}.sorted.bam |bedtools intersect -v -a - -b ${BLACKLIST} > ${sample}.${SZ}.sorted.bed 
rpm=$(wc -l ${sample}.$SZ.sorted.bed |cut -f 1 -d ' '|awk '{printf "%f\n",$1/1000000.0}')
rpm=$(awk '{printf "%f\n",1/("'$rpm'")}' /sibcb1/wuweilab1/wuwei2022/Wuwei/Pipeline/tst.txt)


bedtools genomecov -bg -strand "+" -scale $rpm -i ${sample}.$SZ.sorted.bed -g $CHROMSIZE|bedtools sort -i stdin> ${sample}.${SZ}.pos.bedgraph
bedtools genomecov -bg -strand "-" -scale $rpm -i ${sample}.$SZ.sorted.bed -g $CHROMSIZE|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin > ${sample}.${SZ}.neg.bedgraph
bedGraphToBigWig ${sample}.${SZ}.pos.bedgraph $CHROMSIZE ${sample}.${SZ}.pos.bw
bedGraphToBigWig ${sample}.${SZ}.neg.bedgraph $CHROMSIZE ${sample}.${SZ}.neg.bw
rm ${sample}.${SZ}.pos.bedgraph ${sample}.${SZ}.neg.bedgraph

bedtools genomecov -bg -scale $rpm -i ${sample}.$SZ.sorted.bed -g $CHROMSIZE|bedtools sort -i stdin> ${sample}.${SZ}.bedgraph
bedGraphToBigWig ${sample}.${SZ}.bedgraph $CHROMSIZE ${sample}.${SZ}.bw
rm ${sample}.${SZ}.bedgraph


