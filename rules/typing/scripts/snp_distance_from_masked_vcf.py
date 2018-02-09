import pandas
import itertools
import pysam
import re
import numpy

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

def get_number_of_snps_and_position(snp_list, i, j):
    vect1 = list(snp_list[[i]].values)
    vect2 = list(snp_list[[j]].values)
    diff = 0
    diff_pos = []
    unknowns = 0
    unknowns_pos = []
    for o in range(len(vect1)):
        if snp_list.loc[o, i] == snp_list.loc[o, j]:
            pass
        elif len(snp_list.loc[o, "ALT"].split(","))>1:
            muts=snp_list.loc[o, "ALT"].split(",")
            if snp_list.loc[o, i] and snp_list.loc[o, j]:
                if muts[snp_list.loc[o, i]-1] =="*" and muts[snp_list.loc[o, j]-1] =="*":
                    pass
                elif muts[snp_list.loc[o, i]-1] =="*" or muts[snp_list.loc[o, j]-1] =="*":
                    unknowns += 1
                    unknowns_pos.append(snp_list.loc[o, "POS"])
                elif snp_list.loc[o, i] != snp_list.loc[o, j]:
                    diff += 1
                    diff_pos.append(snp_list.loc[o, "POS"])
            elif snp_list.loc[o, i] and muts[snp_list.loc[o, i]-1]!="*":
                diff += 1
                diff_pos.append(snp_list.loc[o, "POS"])
            elif snp_list.loc[o, j] and muts[snp_list.loc[o, j]-1]!="*":
                diff += 1
                diff_pos.append(snp_list.loc[o, "POS"])
        else:
            if snp_list.loc[o, i] != snp_list.loc[o, j]:
                diff += 1
                diff_pos.append(snp_list.loc[o, "POS"])
    return(diff, diff_pos, unknowns, unknowns_pos)

def get_acgt_count(alignment_file, position):
    samfile = pysam.AlignmentFile(alignment_file, "rb")
    iterator=samfile.pileup(contig=acc, start=position-1, stop=position)
    counts = defaultdict(int)
    for x in iterator:
        if x.reference_pos == (position-1):
            for y in x.pileups:
                if y.query_position is not None:
                    counts[y.alignment.query_sequence[y.query_position]] += 1
    return(counts)

def highlight_max(row):
    values = [ x for x in row if isinstance(x, float) and x<=1 ]
    is_max = values.index(max(values))
    ret = [""]*7
    ret[is_max+2] = 'background-color: yellow' 
    return(ret)

def get_mapping_result_at_position(bam_file, position, sample):
    counts = get_acgt_count(bam_file, position)
    if sum(counts.values()):
        results = [ counts[x]/sum(counts.values()) for x in ["A", "C", "G", "T"] ]
    else:
        results = [0]*4
    return [sample, position] + results + [sum(counts.values())]
    

acc=list(SeqIO.parse(snakemake.input["gbk"], "genbank"))[0].id

snps = pandas.read_csv(snakemake.input["genotype"], sep="\t", header=0)
snps.rename(axis="columns", mapper= lambda x: re.sub("\[[0-9]+\]", "", x.replace(":GT", "").replace("# ", "")), inplace=True)
snps.replace(to_replace=".", value=0, inplace=True)
numeric = snps.drop(["CHROM", "REF", "ALT"], axis=1).apply(pandas.to_numeric).reset_index(drop=True)

pandas.concat([snps.loc[:,["CHROM"]], numeric.loc[:,["POS"]]], axis=1)
snps=pandas.concat([snps.loc[:,["CHROM"]], numeric.loc[:,["POS"]], snps.loc[:, ["REF", "ALT"]].reset_index(drop=True), numeric.drop(["POS"], axis=1)], axis=1)

ref = snakemake.wildcards["ref"]
all_samples = snakemake.params["samples"]
all_samples.sort(key=natural_keys)

number_of_snps = {}
position_of_snps = {}
number_of_unknowns = {}
position_of_unknowns = {}

matrix_distances = pandas.DataFrame(0, index=all_samples + [ref], columns = all_samples + [ref])
matrix_unknowns = pandas.DataFrame(0, index=all_samples + [ref], columns = all_samples + [ref])

for sample1 in all_samples:
    pos = [x for x in snps[[sample1]].values if x != 0]
    diff = len(pos)
    for i in range(len(snps[[sample1]].values)):
        if len(snps.ix[i,"ALT"].split(","))>1 and snps.ix[i, sample1]:
            muts=snps.ix[i,"ALT"].split(",")
            if muts[int(snps.ix[i, sample1])-1]=="*":
                diff -= 1
                matrix_unknowns.loc[sample1, ref] += 1
                matrix_unknowns.loc[ref, sample1] += 1 
    number_of_snps[frozenset((sample1,ref))]=diff
    matrix_distances.loc[sample1, ref] = number_of_snps[frozenset((sample1, ref))]
    matrix_distances.loc[ref, sample1] = number_of_snps[frozenset((sample1, ref))]

for sample1, sample2 in itertools.combinations(all_samples, 2):
    number_of_snps[frozenset((sample1, sample2))], position_of_snps[frozenset((sample1,sample2))], number_of_unknowns[frozenset((sample1, sample2))], position_of_unknowns[frozenset((sample1,sample2))]  = get_number_of_snps_and_position(snps, sample1, sample2)
    matrix_distances.loc[sample1, sample2] = number_of_snps[frozenset((sample1,sample2))]
    matrix_distances.loc[sample2, sample1] = number_of_snps[frozenset((sample1,sample2))]
    matrix_unknowns.loc[sample1, sample2] = number_of_unknowns[frozenset((sample1, sample2))]
    matrix_unknowns.loc[sample2, sample1] = number_of_unknowns[frozenset((sample1, sample2))]
    
bam_files=str(snakemake.input["bams"]).split(" ")

mapping_at_snps = []
for i, j in itertools.combinations(all_samples, 2):
    for k in position_of_snps[frozenset((i,j))]:
        if number_of_snps[frozenset((i,j))] < snakemake.params["dist_thre"]:
            mapping_at_snps.append(get_mapping_result_at_position([s for s in bam_files if i in s][0], k, i))
            mapping_at_snps.append(get_mapping_result_at_position([s for s in bam_files if j in s][0], k, j))

df = pandas.DataFrame(mapping_at_snps, columns=['Strain', 'Position in reference genome', 'A', 'C', 'G', 'T', 'Total Coverage'])

out_csv_distances = snakemake.output["out_csv_dist"]
out_csv_positions = snakemake.output["out_csv_pos"]
out_xlsx_distances = snakemake.output["out_xlsx_dist"]
out_xlsx_positions = snakemake.output["out_xlsx_pos"]

unknowns_number_xlsx = snakemake.output["out_xlsx_unknowns"]
unknowns_number_csv = snakemake.output["out_csv_unknowns"]

matrix_distances.to_csv(out_csv_distances, index= True, sep="\t")

writer = pandas.ExcelWriter(out_xlsx_distances)
matrix_distances.to_excel(writer, snakemake.wildcards["ref"], index=True)
writer.save()

df.to_csv(out_csv_positions, index=False, sep="\t")
    
writer = pandas.ExcelWriter(out_xlsx_positions)
df.style.apply(highlight_max, axis=1).to_excel(writer, ref, index=False, engine="openpyxl")
writer.save()



matrix_unknowns.to_csv(unknowns_number_csv, index=True, sep="\t")
writer = pandas.ExcelWriter(unknowns_number_xlsx)
matrix_unknowns.to_excel(writer, snakemake.wildcards["ref"], index=True)
writer.save()
