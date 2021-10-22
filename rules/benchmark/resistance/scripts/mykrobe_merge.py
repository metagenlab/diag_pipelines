import json
import pandas
import os

mykrobe_results = snakemake.input["all_mykrobe"]

sample2drug2susceptibility = {}

nr_drug_list = []

for mykrobe_result in mykrobe_results:
    sample = mykrobe_result.split("/")[1]
    sample2drug2susceptibility[sample] = {}
    result_table = pandas.read_csv(mykrobe_result, sep=",", index_col="drug")
    for drug in result_table.index:
        drug_l = drug.lower()
        if drug_l not in nr_drug_list:
            nr_drug_list.append(drug_l)
        sample2drug2susceptibility[sample][drug_l] = result_table.at[drug, "susceptibility"]

header = ["sample"] + nr_drug_list
with open(snakemake.output[0], 'w') as f:
    f.write("\t".join(header)+'\n')
    for sample in sample2drug2susceptibility:
        f.write(f"{sample}\t" + '\t'.join([sample2drug2susceptibility[sample][i] for i in nr_drug_list])+'\n')
