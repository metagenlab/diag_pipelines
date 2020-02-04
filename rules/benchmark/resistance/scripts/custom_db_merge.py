#!/usr/bin/env python

import vcf

custom_db_results = snakemake.input["all_custom_db"]
resistance_genes_table = snakemake.input["resistance_genes"]

nr_drug_list = []

gene2resistance = {}
with open(resistance_genes_table, "r") as f:
    for n,row in enumerate(f):
        if n == 0:
            continue
        else:
            data = row.rstrip().split("\t")
            drug = data[0].lower()
            gene = data[1]
            if drug not in nr_drug_list:
                nr_drug_list.append(drug)
            if gene not in gene2resistance:
                gene2resistance[gene] = [drug]
            else:
                gene2resistance[gene].append(drug)



sample2drug2susceptibility = {}

sample2drug2resistance = {}

for custom_db_result in custom_db_results:
    sample = custom_db_result.split("/")[1]
    sample2drug2susceptibility[sample] = {}
    
    vcf_reader = vcf.Reader(open(custom_db_result, 'r'))
    for record in vcf_reader:
        gene = record.INFO["GENE"]
        drugs = gene2resistance[gene]
        for drug in drugs:
            drug_l = drug.lower()
            sample2drug2susceptibility[sample][drug_l] = "R"
    for drug in nr_drug_list:
        if drug not in sample2drug2susceptibility[sample]:
            sample2drug2susceptibility[sample][drug] = "S"
            

header = ["sample"] + nr_drug_list
with open(snakemake.output[0], 'w') as f:

    f.write("\t".join(nr_drug_list)+'\n')
    for sample in sample2drug2susceptibility:
        f.write(f"{sample}\t" + '\t'.join([sample2drug2susceptibility[sample][i] for i in nr_drug_list])+'\n')
