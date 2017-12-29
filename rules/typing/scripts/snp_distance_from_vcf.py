import pandas
import itertools

snps = pandas.read_csv(snakemake.input[0], sep="\t", header=0)

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


with open(snakemake.output[1], "w") as posfile:
    for i, j in itertools.combinations_with_replacement(snakemake.params["samples"], 2):
        if i != j:
            for k in position_of_snps[frozenset((i,j))]:
                posfile.write(str(i)+"-"+str(j)+"\t"+str(k)+"\n")
