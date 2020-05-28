import json
import pandas
import os
import re

rgi_antibiotics = snakemake.input["all_rgi"]

sample2drug2susceptibility = {}

nr_drug_list = []


#ARO	Name	Software	Gene	Resistance Type	Variant	Antibiotic resistance prediction	ARO_antibiotic	Class	Mechanism	AMR family
#ARO:3004628	AAC(2')-IIa	rgi	AAC(2')-IIa	homology model		aminoglycoside	ARO:0000016	aminoglycoside 	antibiotic inactivation	AAC(2')


for rgi_result in rgi_antibiotics:
    sample = rgi_result.split("/")[1]
    sample2drug2susceptibility[sample] = {}
    result_table = pandas.read_csv(rgi_result, sep="\t", index_col="ARO")
    for n, row in result_table.iterrows():
        mechanism = row["Resistance Mechanism"] 
        drug_list = [re.sub(" antibiotic","", i) for i in row["Drug Class"].lower().split("; ")]
        for drug in drug_list:
            if drug not in nr_drug_list:
                nr_drug_list.append(drug)
            if "antibiotic efflux" in mechanism:
                continue
            sample2drug2susceptibility[sample][drug] = "R"

for sample in sample2drug2susceptibility:
    for drug in nr_drug_list:
        if drug not in sample2drug2susceptibility[sample]:
            sample2drug2susceptibility[sample][drug] = "S"

header = ["sample"] + nr_drug_list
with open(snakemake.output[0], 'w') as f:
    f.write("\t".join(header)+'\n')
    for sample in sample2drug2susceptibility:
        f.write(f"{sample}\t" + '\t'.join([sample2drug2susceptibility[sample][i] for i in nr_drug_list])+'\n')
