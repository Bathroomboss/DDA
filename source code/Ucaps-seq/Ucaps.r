#!/usr/bin/env Rscript

###This script uses the sequences of 5bp upstream and downstream of the putative break sites as input to generate ggseqlog plots

args = (commandArgs(TRUE))
seqname <- paste(args[[1]],".breaksite_5bp.seq.txt", sep="")
logoname <- paste(args[[1]],"breaksite_5bp_seqlogo.pdf", sep="",collapse="")

library(ggseqlogo)
library(viridis)
library(Biostrings)
library(ggplot2)

seq<-read.table(seqname)

sequences <- DNAStringSet(seq[,1])
pdf(file=logoname)
ggseqlogo(consensusMatrix(sequences),seq_type='dna')+
theme_bw() +
theme(panel.grid=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.position="none",
        panel.background = element_rect(fill = "#F5F8F9"),
        plot.background = element_rect(fill = "#F5F8F9")) +
        scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11), labels = c("-5", "-4", "-3", "-2", "-1", "Ucaps", "1", "2", "3", "4", "5"))
dev.off()

