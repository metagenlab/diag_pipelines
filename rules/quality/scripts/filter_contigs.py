#!/usr/bin/env python

from Bio import SeqIO
import pandas
import numpy as np

print("OK")


records = SeqIO.parse(snakemake.input[0], "fasta")

header = ["Name", "Length", "Mapped_bases", "Mean_depth", "std"]
df_depth = pandas.read_csv(snakemake.input[1], sep="\t", names=header).set_index("Name")

depth_cutoff = float(snakemake.params["cutoff"])


def N50(list_of_lengths):
    total = sum(list_of_lengths)
    half = float(total) / 2
    cumul = 0
    list_of_lengths.sort(reverse=True)
    for c in list_of_lengths:
        cumul += c
        if cumul >= half:
            return c


n50 = N50(df_depth["Length"].to_list())
largest_contigs = df_depth[df_depth["Length"] >= n50]

average_depth = np.average(largest_contigs["Mean_depth"], 
                           weights=largest_contigs["Length"])

fract_cutoff = average_depth * depth_cutoff

high_depth = []
low_depth = []
for record in records:
    print(record)
    if df_depth.loc[record.name, "Mean_depth"] < fract_cutoff:
        low_depth.append(record)
    else:
        high_depth.append(record)



SeqIO.write(high_depth, snakemake.output[0], "fasta")
SeqIO.write(low_depth, snakemake.output[1], "fasta")
