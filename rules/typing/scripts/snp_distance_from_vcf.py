import pandas
import itertools
import pysam
import re
from collections import defaultdict

from Bio import SeqIO


acc=re.sub("\..*", "", list(SeqIO.parse(snakemake.input["gbk"], "genbank"))[0].id)


snps = pandas.read_csv(snakemake.input["genotype"], sep="\t", header=0)

ref = snakemake.wildcards["ref"]
all_samples = sorted(snakemake.params["samples"] + [ref])
number_of_snps = {}

position_of_snps={}

for i in snakemake.params["samples"]:
    vect = [ 1 if x=="1" else 0 for x in list(snps[[i]].values) ]
    number_of_snps[frozenset((i,ref))]=sum(vect)
    vect = "".join([str(i) for i in vect])
    for j in all_samples:
        if j != i and frozenset((i, j)) not in number_of_snps.keys():
            vect2 = "".join([ "1" if x=="1" else "0" for x in snps[[j]].values ])
            diff = 0
            for x, y in zip(vect, vect2):
                if x != y:
                    diff += 1
            number_of_snps[frozenset((i,j))] = diff
            position_of_snps[frozenset((i,j))] = list(snps.loc[[ i for i in range(len(vect)) if vect[i] != vect2[i] ], "POS"].values)

res_str="\t"+"\t".join(all_samples)+"\n"
for j in all_samples:
    res_str+=j+"\t"
    for i in all_samples:
        if i != j:
            res_str+=str(number_of_snps[frozenset((i,j))])+"\t"
        else:
            res_str+="0\t"
    res_str+="\n"

with open(snakemake.output[0], "w") as resfile:
    resfile.write(res_str)

bam_files=str(snakemake.input["bams"]).split(" ")


def return_acgt_count(alignment_file, position):
    samfile = pysam.AlignmentFile(alignment_file, "rb")
    iterator=samfile.pileup(contig=acc, start=position-1, stop=position)
    counts={}
    for x in ["A", "C", "G", "T"]:
        counts[x]=0
    for x in iterator:
        if x.reference_pos == (position-1):
            for y in x.pileups:
                if y.query_position is not None:
                    counts[y.alignment.query_sequence[y.query_position]] += 1
    return(counts)
        
    


with open(snakemake.output[1], "w") as posfile:
    for i, j in itertools.combinations_with_replacement(snakemake.params["samples"], 2):
        if i != j:
            if number_of_snps[frozenset((i,j))] < snakemake.params["dist_thre"]:
                for k in position_of_snps[frozenset((i,j))]:
                    file1 = [s for s in bam_files if i in s][0]
                    file2 = [s for s in bam_files if j in s][0]
                    c1 = return_acgt_count(file1, k)
                    c2 = return_acgt_count(file2, k)
                    compa = ""
                    for nucl in sorted(c1.keys()):
                        compa += str(c1[nucl])+"-"+str(c2[nucl])+"\t"
                    posfile.write(str(i)+"\t"+str(j)+"\t"+str(k)+"\t"+compa+"\n")
