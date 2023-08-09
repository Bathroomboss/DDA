#!/bin/bash
# HTGTS

##############################################################################
# Note:
# If there are result files on GEO, we directly download and extract prey information to generate a BED file. 
# Otherwise, we primarily follow the instructions on this website for analysis: https://robinmeyers.github.io/transloc_pipeline/thedocs.html
# Generate a result.tlx file and then extract prey information to generate a BED file.
##############################################################################
# example:
# go to GEO database and find Series record "ftp"

cutadapt  -e 0.1 -O 3 -m 55 --quality-cutoff 25 -a AGATCGGAAGAG -A AGATCGGAAGAG \ 
-o cutadapt/${Filename}_cutadapt.R1_fastq.gz -p cutadapt/${Filename}_cutadapt_R2.fastq.gz  ${Filename}_1.fq.gz ${Filename}_2.fq.gz
TranslocWrapper.pl metadata.txt cutadapt/ results/ --threads 8  --pipeline-opt --no-dedup 

awk 'BEGIN{FS=OFS="\t"}{if ($5==1) print $3,$6,$7,$1,$18,"+";else print $3,$6,$7,$1,$18,"-"}'  $File.txt > $File.bed

sed -i '1d' $File.bed
sort -k1,1 -k2,2n $File.bed |bedtools intersect -a - -b $Blacklist -v > $File.sorted.bed
bedtools genomecov -bg -strand "+" -i ${File}.sorted.bed -g ${chromSize}|bedtools sort -i stdin> ${File}.pos.bedgraph
bedtools genomecov -bg -strand "-" -i ${File}.sorted.bed -g ${chromSize}|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin > ${File}.neg.bedgraph
bedGraphToBigWig ${File}.pos.bedgraph ${chromSize} ${File}.pos.bw
bedGraphToBigWig ${File}.neg.bedgraph ${chromSize} ${File}.neg.bw
rm ${File}.pos.bedgraph ${File}.neg.bedgraph 