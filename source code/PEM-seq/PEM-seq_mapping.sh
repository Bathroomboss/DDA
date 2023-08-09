#!/bin/bash
# PEM-seq

##############################################################################
# Note:
# If there are result files on GEO, we directly download and extract prey information to generate a BED file. 
# Otherwise, we primarily follow the instructions on this website for analysis: https://github.com/liumz93/PEM-Q
# Generate a Translocation.tab file and then extract prey information to generate a BED file.
##############################################################################

PEM-Q.py <genome> <sample> <cutsite> <primer_chr> <primer_start> <primer_end> <primer_strand> <primer>

awk 'BEGIN{FS=OFS="\t"}{if ($7=="+") print $6,$8,$9,$1,".","+";else print $6,$8,$9,$1,".","-"}'  ${sample}_Translocation.tab > ${sample}.bed

sed -i '1d' ${sample}.bed
sort -k1,1V -k2,2n -k3,3n ${sample}.bed |bedtools intersect -a - -b $BLACKLIST -v > ${sample}.${SZ}.sorted.bed
bedtools genomecov -bg -strand "+" -i ${sample}.${SZ}.sorted.bed -g ${CHROMSIZE}|bedtools sort -i stdin > ${sample}.${SZ}.pos.bedgraph
bedtools genomecov -bg -strand "-" -i ${sample}.${SZ}.sorted.bed -g ${CHROMSIZE}|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin > ${sample}.${SZ}.neg.bedgraph
bedtools genomecov -bg -i ${sample}.${SZ}.sorted.bed -g ${CHROMSIZE}|bedtools sort -i stdin > ${sample}.${SZ}.bedgraph
bedGraphToBigWig ${sample}.${SZ}.pos.bedgraph ${CHROMSIZE} ${sample}.${SZ}.pos.bw
bedGraphToBigWig ${sample}.${SZ}.neg.bedgraph ${CHROMSIZE} ${sample}.${SZ}.neg.bw
bedGraphToBigWig ${sample}.${SZ}.bedgraph ${CHROMSIZE} ${sample}.${SZ}.bw
rm ${sample}.${SZ}.pos.bedgraph  ${sample}.${SZ}.neg.bedgraph ${sample}.${SZ}.bedgraph

