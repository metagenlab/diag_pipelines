import json
import pandas
import os
import re
import math 

mykrobe_results = snakemake.input["all_mykrobe"]

sample2drug2variant = {}

for mykrobe_result in mykrobe_results:
    sample = mykrobe_result.split("/")[1]
    sample2drug2variant[sample] = {}
    result_table = pandas.read_csv(mykrobe_result, sep="\t", index_col="drug")
    for drug in result_table.index:
        drug_l = drug.lower()
        if drug_l not in sample2drug2variant[sample]:
            sample2drug2variant[sample][drug_l] = []
        
        # rrs_A1401X-A1473246G:152:0:16577
        # fabG1_C-15X-C1673425T:180:0:21468
        variants = result_table.at[drug, "variants (gene:alt_depth:wt_depth:conf)"]
        if isinstance(variants, str):
            variants = variants.split(";")
            for variant in variants:
                regex = '([A-Za-z\d]+)_([A-Za-z\d\-]+)-([A-Za-z\*:\d]+)'
                s = re.search(regex, variant)
                gene = s.group(1)
                change = s.group(2)
                variant_format = "%s##%s" % (gene, change)
                sample2drug2variant[sample][drug_l].append(variant_format)
        #genes = result_table.at[drug, "genes (prot_mut-ref_mut:percent_covg:depth)"]


with open(snakemake.output[0], "w") as f:
    for sample in sample2drug2variant:
        for drug in sample2drug2variant[sample]:
            for variant in sample2drug2variant[sample][drug]:
                # rrs_A1401X-A1473246G:152:0:16577
                # fabG1_C-15X-C1673425T:180:0:21468
                gene, change = variant.split("##")
                regex = '([A-Za-z]+)([\-\d]+)([A-Za-z]+)'
                s = re.search(regex, change)
                REF = s.group(1)
                position = s.group(2)
                ALT = s.group(3)

                if len(REF) < len(ALT):
                    vartype = 'INS'
                elif len(REF) > len(ALT):
                    vartype = 'DEL'
                else:
                    vartype = 'SNP'
                
                f.write(f"mykrobe\t{sample}\t{drug}\t{gene}\t{position}\t{change}\t{vartype}\n")
