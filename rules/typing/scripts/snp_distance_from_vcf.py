import pandas
import itertools
import pysam
import re

from Bio import SeqIO

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]

out_xlsx = snakemake.output[2]
out_csv = snakemake.output[0]

acc=re.sub("\..*", "", list(SeqIO.parse(snakemake.input["gbk"], "genbank"))[0].id)

snps = pandas.read_csv(snakemake.input["genotype"], sep="\t", header=0)

ref = snakemake.wildcards["ref"]
all_samples = snakemake.params["samples"]
all_samples.sort(key=natural_keys)
number_of_snps = {}
position_of_snps = {}

result = pandas.DataFrame(0, index=all_samples + [ref], columns = all_samples + [ref])

for i, j in itertools.combinations(all_samples, 2):
    vect = list(snps[[i]].values)
    if frozenset((i, ref)) not in number_of_snps.keys():
        number_of_snps[frozenset((i, ref))] = sum(vect)[0]
        result.loc[[i], [ref]] = number_of_snps[frozenset((i, ref))]
        result.loc[[ref], [i]] = number_of_snps[frozenset((i, ref))]
    vect2 = list(snps[[j]].values)
    if frozenset((j, ref)) not in number_of_snps.keys():
        number_of_snps[frozenset((j, ref))] = sum(vect2)[0]
        result.loc[[j], [ref]] = number_of_snps[frozenset((j, ref))]
        result.loc[[ref], [j]] = number_of_snps[frozenset((j, ref))]
    number_of_snps[frozenset((i,j))] = len([u for u in range(len(vect)) if vect[u] != vect2[u]])
    position_of_snps[frozenset((i,j))] = list(snps.loc[[ u for u in range(len(vect)) if vect[u] != vect2[u] ], "POS"].values)
    result.loc[[i], [j]] = number_of_snps[frozenset((i,j))]
    result.loc[[j], [i]] = number_of_snps[frozenset((i,j))]
    
writer = pandas.ExcelWriter(out_xlsx)
result.to_excel(writer, snakemake.wildcards["ref"], index=True)
writer.save()

result.to_csv(out_csv, index= True, sep="\t")
writer.save()
    
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
