import json
import pandas
import os
import re

rgi_antibiotics = snakemake.input["all_rgi"]

nr_drug_list = []


#ARO	Name	Software	Gene	Resistance Type	Variant	Antibiotic resistance prediction	ARO_antibiotic	Class	Mechanism	AMR family
#ARO:3004628	AAC(2')-IIa	rgi	AAC(2')-IIa	homology model		aminoglycoside	ARO:0000016	aminoglycoside 	antibiotic inactivation	AAC(2')

sample2drug2variant = {}

for rgi_result in rgi_antibiotics:
    result_table = pandas.read_csv(rgi_result, sep="\t", index_col="ARO")
    sample = rgi_result.split("/")[1]
    sample2drug2variant[sample] = {}
    for n, row in result_table.iterrows():
        mechanism = row["Resistance Mechanism"].split(", ")
        model_type = row["Model_type"]
        snps = row["SNPs_in_Best_Hit_ARO"]
        
        aro_name = row["Best_Hit_ARO"]
        drug_list = [re.sub(" antibiotic","", i) for i in row["Drug Class"].lower().split("; ")]
        
        # resistance due to SNPs
        if isinstance(snps, str):
            snps = snps.split(", ")
            for drug in drug_list:
                if "antibiotic efflux" in mechanism:
                    continue
                if drug not in sample2drug2variant[sample]:
                    sample2drug2variant[sample][drug] = []
                for snp in snps:
                    var_label = "%s##%s##%s" % (aro_name, model_type, snp)
                    sample2drug2variant[sample][drug].append(var_label)
        
        # resistance NOT due to SNPs
        else:
            for drug in drug_list:
                if "antibiotic efflux" in mechanism:
                    continue
                if drug not in sample2drug2variant[sample]:
                    sample2drug2variant[sample][drug] = []

                var_label = "%s##%s##%s" % (aro_name, model_type, "nosnp")
                sample2drug2variant[sample][drug].append(var_label)

with open(snakemake.output[0], 'w') as f:
    for sample in sample2drug2variant:
        for drug in sample2drug2variant[sample]:
            for variant in sample2drug2variant[sample][drug]:
                gene, model, change = variant.split("##")
                regex = '[A-Za-z]+([\-\d]+)[A-Za-z]+'
                if change != "nosnp":
                    s = re.search(regex, change)
                    position = s.group(1)
                    vartype = "SNP"
                else:
                    position = None
                    vartype = "nosnp"
                f.write(f"rgi\t{sample}\t{drug}\t{gene}\t{position}\t{change}\t{vartype}\n")
