import pandas
import csv
import re
import itertools

def get_number_of_snps_and_position(snp_list, i, j):
    vect1 = list(snp_list[[i]].values)
    vect2 = list(snp_list[[j]].values)
    return( len([u for u in range(len(vect1)) if vect1[u] != vect2[u]]), list(snp_list.loc[[ u for u in range(len(vect1)) if vect1[u] != vect2[u] ], "POS"].values))

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]

all_samples=snakemake.params["samples"]
all_samples.sort(key=natural_keys)
ref = snakemake.wildcards["ref"]

snps = pandas.read_csv(snakemake.input["genotype"], sep="\t", header=0)
snps.rename(axis="columns", mapper= lambda x: re.sub("\[[0-9]+\]", "", x.replace(":GT", "")), inplace=True)
snps.replace(to_replace=".", value=0, inplace=True)

evidence_files=str(snakemake.input["evidences"]).split(" ")

number_of_snps = {}


files = {}
for i in all_samples:
    files[i] = pandas.read_csv([s for s in evidence_files if i in s][0], sep="\t", header=None, names=["POS", "REF", "COV", "A", "C", "G", "T"], usecols=[1,2,3,4,5,6,7], index_col=0)
    files[i] = files[i].transform(lambda x: [int(u.split(":")[1]) if isinstance(u, str) and len(u.split(":"))>1 else u for u in x])
    files[i]["MAX"] = files[i][["A", "C", "G", "T"]].idxmax(axis=1)

for sample1 in all_samples:
    diff = len([ x for x in snps[[sample1]].values if x != 0])
    number_of_snps[frozenset((sample1,ref))]=diff

    
for sample1, sample2 in itertools.combinations(all_samples, 2):
    number_of_snps[frozenset((sample1, sample2))], position_of_snps = get_number_of_snps_and_position(snps, sample1, sample2)
    for k in position_of_snps:
        if k not in files[sample1].index or k not in files[sample2].index:
            number_of_snps[frozenset((sample1, sample2))] -= 1
        elif files[sample1].loc[k,"COV"] < 10 or files[sample2].loc[k,"COV"] < 10:
            number_of_snps[frozenset((sample1, sample2))] -= 1
        elif files[sample1].loc[k, "MAX"] == files[sample2].loc[k, "MAX"]:
            if files[sample1].loc[k, files[sample1].loc[k, "MAX"]]/files[sample1].loc[k, "COV"] > 0.8 and files[sample2].loc[k, files[sample1].loc[k, "MAX"]]/files[sample2].loc[k, "COV"] > 0.8:
                number_of_snps[frozenset((sample1, sample2))] -= 1

            




matrix_distances = pandas.DataFrame(0, index=all_samples + [ref], columns = all_samples + [ref])

for sample1 in all_samples:
    matrix_distances.loc[sample1, ref] = number_of_snps[frozenset((sample1, ref))]
    matrix_distances.loc[ref, sample1] = number_of_snps[frozenset((sample1, ref))]
    
for sample1, sample2 in itertools.combinations(all_samples, 2):
    matrix_distances.loc[sample1, sample2] = number_of_snps[frozenset((sample1,sample2))]
    matrix_distances.loc[sample2, sample1] = number_of_snps[frozenset((sample1,sample2))]


out_xlsx_distances = snakemake.output["out_xlsx_dist"]
out_csv_distances = snakemake.output["out_csv_dist"]



matrix_distances.to_csv(out_csv_distances, index= True, sep="\t")

writer = pandas.ExcelWriter(out_xlsx_distances)
matrix_distances.to_excel(writer, snakemake.wildcards["ref"], index=True)
writer.save()



