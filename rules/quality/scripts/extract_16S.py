
from Bio import SeqIO
records = SeqIO.parse(snakemake.input[0], "fasta") 
keep = []
for record in records:
    if "16S_rRNA" in record.description:
        keep.append(record)
SeqIO.write(keep, snakemake.output[0], "fasta")