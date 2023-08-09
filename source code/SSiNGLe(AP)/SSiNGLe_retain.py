#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
from io import StringIO
def read_fasta_file(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
        seqs = []
        seq = None
        for line in lines:
            if line.startswith('>'):
                if seq is not None:
                    seqs.append(seq)
                seq = [line.strip(), '']
            else:
                seq[1] += line.strip()
        if seq is not None:
            seqs.append(seq)
    return seqs

def filter_seqs(seqs, min_a_count):
    new_seqs = []
    for s in seqs:
        seq_title = s[0].replace(">", '').replace(":", "\t").replace("-", "\t")
        seq = s[1].lower()
        a_count = seq.count('a')
        if a_count <=  min_a_count:
            new_seqs.append([seq_title, seq])
    return new_seqs

def write_fasta_file(seqs, output_filename):
    with open(output_filename, 'w') as f:
        for s in seqs:
            f.write('{}\n{}\n'.format(s[0], s[1]))

def main():
    parser = argparse.ArgumentParser(description='screene A enrichment sequence')
    parser.add_argument('input_filename', type=str, help='input fasta file')
    parser.add_argument('output_filename', type=str, help='output fasta file')
    parser.add_argument('-c', '--min-a-count', type=int, default=8, help='Minimum A counts,defult=8')

    args = parser.parse_args()

    seqs = read_fasta_file(args.input_filename)
    seqs = filter_seqs(seqs, args.min_a_count)
    write_fasta_file(seqs, args.output_filename)


if __name__ == '__main__':
    main()
