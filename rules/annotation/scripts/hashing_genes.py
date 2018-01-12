import hashlib
from Bio import SeqIO

with open(snakemake.output[0], "w") as hash_file:
    for record in SeqIO.parse(snakemake.input[0], "fasta"):
        hash_file.write(hashlib.sha256(str(record.seq).encode()).hexdigest()+"\n")
