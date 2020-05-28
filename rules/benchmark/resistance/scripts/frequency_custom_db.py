#!/usr/bin/env python

import vcf
import re 

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

sample2drug2variant = {}


for custom_db_result in custom_db_results:
    sample = custom_db_result.split("/")[1]
    sample2drug2variant[sample] = {}
    vcf_reader = vcf.Reader(open(custom_db_result, 'r'))
    for record in vcf_reader:
        gene = record.INFO["GENE"]
        wildtype_aa = record.INFO["WILDTYPE_AA_NUCL"]
        mutated_aa = record.INFO["MUTATED_AA_NUCL"]
        position = record.INFO["POSITION_IN_GENE"]
        drugs = gene2resistance[gene]
        label = "%s##%s%s%s" % (gene, wildtype_aa, position, mutated_aa)
        for drug in drugs:
            drug_l = drug.lower()
            if drug_l not in sample2drug2variant[sample]:
                sample2drug2variant[sample][drug_l] = [label]
            else:
                sample2drug2variant[sample][drug_l].append(label)
            
with open(snakemake.output[0], 'w') as f:
    for sample in sample2drug2variant:
        for drug in sample2drug2variant[sample]:
            for variant in sample2drug2variant[sample][drug]:
                gene, change = variant.split("##")
                regex = '[A-Za-z]+([\-\d]+)[A-Za-z\*]+'
                s = re.search(regex, change)
                position = s.group(1)
                f.write(f"custom\t{sample}\t{drug}\t{gene}\t{position}\t{change}\tSNP\n")
