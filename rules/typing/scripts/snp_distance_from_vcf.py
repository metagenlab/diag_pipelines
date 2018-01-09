import pandas
import itertools
import pysam
import re

from collections import defaultdict
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
        result.loc[i, ref] = number_of_snps[frozenset((i, ref))]
        result.loc[ref, i] = number_of_snps[frozenset((i, ref))]
    vect2 = list(snps[[j]].values)
    if frozenset((j, ref)) not in number_of_snps.keys():
        number_of_snps[frozenset((j, ref))] = sum(vect2)[0]
        result.loc[j, ref] = number_of_snps[frozenset((j, ref))]
        result.loc[ref, j] = number_of_snps[frozenset((j, ref))]
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
    counts = defaultdict(int)
    for x in iterator:
        if x.reference_pos == (position-1):
            for y in x.pileups:
                if y.query_position is not None:
                    counts[y.alignment.query_sequence[y.query_position]] += 1
    return(counts)




res = []
with open(snakemake.output[1], "w") as posfile:
    for i, j in itertools.combinations(all_samples, 2):
        if number_of_snps[frozenset((i,j))] < snakemake.params["dist_thre"]:
            for k in position_of_snps[frozenset((i,j))]:
                file1 = [s for s in bam_files if i in s][0]
                file2 = [s for s in bam_files if j in s][0]
                c1 = return_acgt_count(file1, k)
                c2 = return_acgt_count(file2, k)
                res1 = [c1[x]/sum(c1.values()) for x in ["A", "C", "G", "T"]]
                res2 = [c2[x]/sum(c2.values()) for x in ["A", "C", "G", "T"]]
                res.append([i, k, res1[0], res1[1], res1[2], res1[3], sum(c1.values())])
                res.append([j, k, res2[0], res2[1], res2[2], res2[3], sum(c2.values())])



def color_negative_red(val):
    """
    Takes a scalar and returns a string with
    the css property `'color: red'` for negative
    strings, black otherwise.
    """
    if isinstance(val, str) or val < 0.5 :
        color = 'black'
    else:
        color = 'red'
    return 'color: %s' % color


print(pandas.__file__)

df = pandas.DataFrame(res, columns=['Strain', 'Position in reference genome', 'A', 'C', 'G', 'T', 'Total Coverage'])

out_xlsx_positions = snakemake.output[3]
out_csv_positions = snakemake.output[1]

writer = pandas.ExcelWriter(out_xlsx_positions)

df.style.applymap(color_negative_red).to_excel(writer, ref, index=False, engine="openpyxl")
writer.save()

df.to_csv(out_csv_positions, index=False, sep="\t")
