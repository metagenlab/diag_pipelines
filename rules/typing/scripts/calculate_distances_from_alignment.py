from Bio import AlignIO
import itertools


align = AlignIO.read(snakemake.input[0], "fasta")

with open(snakemake.output[0], "w") as dist_file:
    for (i, j) in itertools.combinations(range(len(align)), 2):
        mutations = 0
        pos=0
        print(align[i].id, align[j].id)
        for seq1,seq2 in zip(align[i].seq, align[j].seq):
            pos += 1
            if seq1 != seq2 and seq1 !="N" and seq2 !="N":
                mutations += 1
                if align[i].id=="S7" and align[j].id=="S4":
                    print(pos, mutations, seq1,seq2)
        dist_file.write(align[i].id+" "+align[j].id+" "+ str(mutations)+"\n")
