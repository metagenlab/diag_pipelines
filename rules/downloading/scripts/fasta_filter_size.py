#!/usr/bin/python
from Bio import SeqIO
import sys
import os

#usage: python long.seq.py in.fasta out.fasta 200

fasta_file = snakemake.input[0]
size_limit = snakemake.params[0]
out_file_1 = snakemake.output[0]
out_file_2 = snakemake.output[1]

output_small = open(out_file_1, "w")
output_large = open(out_file_2, "w")

small = []
large = []
for record in SeqIO.parse(open(fasta_file, "rU"), "fasta"):
    if len(record.seq) <= size_limit:
        small.append(record)
    else:
        large.append(record)

SeqIO.write(small, out_file_1, "fasta")
SeqIO.write(large, out_file_2, "fasta")
output_small.close()
output_large.close()