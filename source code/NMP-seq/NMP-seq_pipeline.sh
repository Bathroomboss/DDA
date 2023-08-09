#!/bin/bash
# NMP-seq_mapping
# This pipeline is designed by Qingqing yuan
cat filename.txt | while read line
do
Filename=$line
SZ=sacCer3
File=$Filename.$SZ
Bowtie2Index=/sibcb1/wuweilab1/wuwei-lab/Reference/Yeast/sacCer3/Sequence/Bowtie2Index/genome
Blacklist=/sibcb1/wuweilab1/liangyu/Pipeline/noblcaklist.bed
chromSize=/sibcb1/wuweilab1/wuwei-lab/Reference/Yeast/sacCer3/Sequence/WholeGenomeFasta/genome.chrom.sizes.txt
fa=/sibcb1/wuweilab1/wuwei-lab/Reference/Yeast/sacCer3/Sequence/WholeGenomeFasta/genome.fa

echo $File
bowtie2 -p 24 -x $Bowtie2Index -U $Filename.fastq.gz -S $File.sam 
samtools view -bSh -@ 24 -q 30 $File.sam > $File.bam
samtools sort -@ 24 -o $File.sorted.bam $File.bam
bamToBed -i $File.sorted.bam > $File.sorted.raw.bed 
rm $File.sam $File.bam 

perl -alne 'if($F[5]=~s/\+/-/){$F[2]=$F[1];$F[1]=$F[1]-1;}else{$F[5]=~s/-/\+/;$F[1]=$F[2];$F[2]=$F[2]+1}print join "\t", @F' < $File.sorted.raw.bed > $File.1bp.temp.bed ########

awk '{if ($2 >0) print $0}' $File.1bp.temp.bed > $File.1bp.bed 
rm $File.1bp.temp.bed 

sort -k1,1V -k2,2n -k3,3n $File.1bp.bed >$File.1bp.sorted.bed
rm $File.1bp.bed

bedtools getfasta -s -fi $fa -bed $File.1bp.sorted.bed -fo $File.1bp.strand.fa
awk '{if (NR%2==1) print $0}' $File.1bp.strand.fa > $File.pos.info
awk '{if (NR%2==0) print $0}' $File.1bp.strand.fa > $File.NT.info
python /sibcb2/wuweilab2/wuwei-lab/Data/DamageView/RawData/NMP-seq/NMP_filtering.py $File.pos.info $File.NT.info > $File.filter.temp.bed
rm $File.pos.info $File.NT.info

bedtools intersect -a $File.filter.temp.bed -b $Blacklist -v > $File.filtered.sorted.bed
rm $File.filter.temp.bed

awk '{if ($5=="G"){print $0}}' ${File}.filtered.sorted.bed > ${File}.filtered.sorted_G.bed
awk '{if ($5=="A"){print $0}}' ${File}.filtered.sorted.bed > ${File}.filtered.sorted_A.bed

sort -k1,1V -k2,2n -k3,3n ${File}.filtered.sorted_G.bed >${File}_GREADS.bed
sort -k1,1V -k2,2n -k3,3n ${File}.filtered.sorted_A.bed >${File}_AREADS.bed
rm ${File}.filtered.sorted_G.bed ${File}.filtered.sorted_A.bed

rpm=$(wc -l ${File}_GREADS.bed|cut -f 1 -d ' '|awk '{printf "%f\n",$1/1000000.0}')
rpm=$(awk '{printf "%f\n",1/("'$rpm'")}' /sibcb1/wuweilab1/wuwei2022/Wuwei/Pipeline/tst.txt)

bedtools genomecov -bg -strand "+" -scale $rpm -i ${File}_GREADS.bed -g ${chromSize}|bedtools sort -i stdin> ${File}_GREADS.pos.bedgraph
bedtools genomecov -bg -strand "-" -scale $rpm -i ${File}_GREADS.bed -g ${chromSize}|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin > ${File}_GREADS.neg.bedgraph

bedGraphToBigWig ${File}_GREADS.pos.bedgraph ${chromSize} ${File}_GREADS.pos.bw
bedGraphToBigWig ${File}_GREADS.neg.bedgraph ${chromSize} ${File}_GREADS.neg.bw

rpm=$(wc -l ${File}_AREADS.bed|cut -f 1 -d ' '|awk '{printf "%f\n",$1/1000000.0}')
rpm=$(awk '{printf "%f\n",1/("'$rpm'")}' /sibcb1/wuweilab1/wuwei2022/Wuwei/Pipeline/tst.txt)

bedtools genomecov -bg -strand "+" -scale $rpm -i ${File}_AREADS.bed -g ${chromSize}|bedtools sort -i stdin> ${File}_AREADS.pos.bedgraph
bedtools genomecov -bg -strand "-" -scale $rpm -i ${File}_AREADS.bed -g ${chromSize}|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin > ${File}_AREADS.neg.bedgraph

bedGraphToBigWig ${File}_AREADS.pos.bedgraph ${chromSize} ${File}_AREADS.pos.bw
bedGraphToBigWig ${File}_AREADS.neg.bedgraph ${chromSize} ${File}_AREADS.neg.bw


rm ${File}_GREADS.pos.bedgraph ${File}_GREADS.neg.bedgraph
rm ${File}_AREADS.pos.bedgraph ${File}_AREADS.neg.bedgraph

done
