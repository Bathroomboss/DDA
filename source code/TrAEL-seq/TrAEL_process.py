#!/usr/bin/env python
# -*- coding: utf-8 -*-
#!/usr/bin/python
# -- coding: utf-8 --

import argparse
parser = argparse.ArgumentParser(description='Trim "T" from the beginning of reads in a FASTQ file')
parser.add_argument('infile', help='Input FASTQ file')
parser.add_argument('outfile', help='Output FASTQ file')

args = parser.parse_args()
infile = args.infile
outfile = args.outfile

with open(infile, 'r') as f_in, open(outfile, 'w') as f_out:
    for line in f_in:
        if line.startswith('@'):
            seq_id = line.strip()
            seq = next(f_in).strip()
            if not seq.startswith('T'):
               continue
            original_seq_length = len(seq)
            for i in range(3):
                if seq.startswith('T'):
                    seq = seq[1:]
                else:
                    break
            new_seq_length = len(seq)
            qual_id = next(f_in).strip()
            quality_scores = next(f_in).strip()
            quality_length_difference = original_seq_length - new_seq_length
            updated_qual_id = qual_id.split()[0]
            f_out.write(seq_id + '\n' + seq + '\n' + updated_qual_id + '\n' + quality_scores[quality_length_difference:] + '\n')

