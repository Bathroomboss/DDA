#!/bin/bash
# CPDseq

echo $File
date
if [ ${END} == "SE" ];then
	bowtie2 -p 24 -x $BOWTIE2INDEX -U $DPATH/${File}.fastq.gz -S $DPATH/$File.$SZ.sam --no-unal
else
	bowtie2 -p 24 -x $BOWTIE2INDEX -1 $DPATH/${File}_1.fastq.gz -2 $DPATH/${File}_2.fastq.gz -S $DPATH/$File.$SZ.sam --no-unal
fi

samtools view -bSh -@ 24 -q 30 $DPATH/$File.$SZ.sam > $DPATH/$File.$SZ.bam
rm $DPATH/$File.$SZ.sam
samtools sort -@ 24 -o $DPATH/$File.$SZ.sorted.bam $DPATH/$File.$SZ.bam
rm $DPATH/$File.$SZ.bam
bamToBed -i $DPATH/$File.$SZ.sorted.bam > $DPATH/$File.$SZ.sorted.raw.bed
perl -alne 'if($F[5]=~s/\+/-/){$F[2]=$F[1];$F[1]=$F[1]-2;}else{$F[5]=~s/-/\+/;$F[1]=$F[2];$F[2]=$F[2]+2}print join "\t", @F' < $DPATH/$File.$SZ.sorted.raw.bed > $DPATH/$File.$SZ.2bp.temp.bed
awk '{if ($2 >0) print $0}' $DPATH/$File.$SZ.2bp.temp.bed > $DPATH/$File.$SZ.2bp.bed 
rm $DPATH/$File.$SZ.2bp.temp.bed 
bedtools getfasta -s -fi $FA -bed $DPATH/$File.$SZ.2bp.bed -fo $DPATH/$File.$SZ.2bp.strand.fa
awk '{if (NR%2==1) print $0}' $DPATH/$File.$SZ.2bp.strand.fa > $DPATH/$File.$SZ.pos.info
awk '{if (NR%2==0) print $0}' $DPATH/$File.$SZ.2bp.strand.fa > $DPATH/$File.$SZ.NT.info
python /sibcb2/wuweilab2/liangyu/CPD-seq/CPD_filtering.py $DPATH/$File.$SZ.pos.info $DPATH/$File.$SZ.NT.info > $DPATH/$File.$SZ.filtered.temp.bed
rm $DPATH/$File.$SZ.pos.info $DPATH/$File.$SZ.NT.info
bedtools intersect -a $DPATH/$File.$SZ.filtered.temp.bed -b $BLACKLIST -v > $DPATH/$File.$SZ.filtered.sorted.bed

rm $DPATH/$File.$SZ.filtered.temp.bed
rpm=$(wc -l $DPATH/$File.$SZ.filtered.sorted.bed|cut -f 1 -d ' '|awk '{printf "%f\n",$1/1000000.0}')
rpm=$(awk '{printf "%f\n",1/("'$rpm'")}' /sibcb1/wuweilab1/wuwei2022/Wuwei/Pipeline/tst.txt)
bedtools genomecov -bg -strand "+" -scale $rpm -i $DPATH/$File.$SZ.filtered.sorted.bed -g $CHROMSIZE|bedtools sort -i stdin> $DPATH/$File.$SZ.pos.bedgraph
bedtools genomecov -bg -strand "-" -scale $rpm -i $DPATH/$File.$SZ.filtered.sorted.bed -g $CHROMSIZE|awk '{printf "%s\t%s\t%s\t%s\n",$1,$2,$3,-1*$4}'|bedtools sort -i stdin> $DPATH/$File.$SZ.neg.bedgraph
bedGraphToBigWig $DPATH/$File.$SZ.pos.bedgraph $CHROMSIZE $DPATH/$File.$SZ.pos.bw
bedGraphToBigWig $DPATH/$File.$SZ.neg.bedgraph $CHROMSIZE $DPATH/$File.$SZ.neg.bw
rm $DPATH/$File.$SZ.pos.bedgraph $DPATH/$File.$SZ.neg.bedgraph

date


