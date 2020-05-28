import pandas
import Bio
from Bio.Data.SCOPData import protein_letters_3to1
import re

drug_list =["rifampicin", "isoniazid", "pyrazinamide", 
          "ethambutol", "streptomycin", "fluoroquinolones", "moxifloxacin", 
          "ofloxacin", "levofloxacin", "ciprofloxacin", "aminoglycosides", "amikacin", 
          "kanamycin", "capreomycin", "ethionamide", "para-aminosalicylic_acid", 
          "cycloserine", "linezolid", "bedaquiline", "clofazimine", "delamanid"]  

table_a = pandas.read_csv(snakemake.input[0] , index_col=0, sep="\t", header=0)

sample2drug2variant = {}

for sample, row in table_a.iterrows():
    sample2drug2variant[sample] = {}
    for drug in drug_list:
        sample2drug2variant[sample][drug] = [] 
        if table_a.at[sample, drug] == '-':
            continue
        else:
            variant_list = [i.strip() for i in table_a.at[sample, drug].split(",")]
            for variant in variant_list:
                sample2drug2variant[sample][drug].append(variant)


with open(snakemake.output[0], "w") as f:
    for sample in sample2drug2variant:
        for drug in sample2drug2variant[sample]:
            for variant in sample2drug2variant[sample][drug]:
                # rpoB_p.Leu430Pro
                gene, change = variant.split("_", 1)
                if change[0] == 'p':
                    regex = 'p.([A-Za-z]+)([\d\-]+)([A-Za-z\*]+)'
                    s = re.search(regex, change)
                    REF = protein_letters_3to1[s.group(1).upper()]
                    pos = s.group(2)
                    ALT = s.group(3).upper()
                    if ALT != "*":
                        ALT = protein_letters_3to1[ALT]
                        vartype = "SNP"
                    else:
                        vartype = "STOP"            
                elif change[0] in ['r', 'c']:
                    # r.1401a>g
                    # c.-15C>T
                    # c.470_471insA
                    if 'ins' in change:
                        data = change[2:len(change)].split("ins")
                        pos = data[0]
                        REF = None
                        ALT = data[1]
                        vartype = "INS" 
                    elif 'del' in change:
                        data = change[2:len(change)].split("del")
                        pos = data[0]
                        REF = None
                        ALT = None
                        vartype = "DEL" 
                    else:
                        regex = '[rc].([\d\-]+)([A-Za-z]+)>([A-Za-z]+)'
                        s = re.search(regex, change)
                        pos = s.group(1)
                        REF = s.group(2).upper()
                        ALT = s.group(3).upper()      
                        vartype = "SNP" 
                elif "Chromosome" in change:
                    # gid Chromosome:g.4407986_4408185del
                    if 'ins' in change:
                        print(change)
                        data = change.split("g.")[1].split("ins")
                        pos = data[0]
                        REF = None
                        ALT = data[1]
                        vartype = "chr_INS"
                    elif 'del' in change:
                        print(change)
                        data = change.split("g.")[1].split("del")
                        pos = data[0]
                        REF = None
                        ALT = None
                        vartype = "chr_DEL"
                else:
                    print(gene, change)
                    raise IOError("error")       
                f.write(f"tbprofiler\t{sample}\t{drug}\t{gene}\t{pos}\t{REF}{pos}{ALT}\t{vartype}\n")          

            
