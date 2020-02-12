import pandas

header =["sample", "main_lineage", "sub_lineage", "DR_type", 
          "MDR", "XDR", "rifampicin", "isoniazid", "pyrazinamide", 
          "ethambutol", "streptomycin", "fluoroquinolones", "moxifloxacin", 
          "ofloxacin", "levofloxacin", "ciprofloxacin", "aminoglycosides", "amikacin", 
          "kanamycin", "capreomycin", "ethionamide", "para-aminosalicylic_acid", 
          "cycloserine", "linezolid", "bedaquiline", "clofazimine", "delamanid"] 

drug_list =["rifampicin", "isoniazid", "pyrazinamide", 
          "ethambutol", "streptomycin", "fluoroquinolones", "moxifloxacin", 
          "ofloxacin", "levofloxacin", "ciprofloxacin", "aminoglycosides", "amikacin", 
          "kanamycin", "capreomycin", "ethionamide", "para-aminosalicylic_acid", 
          "cycloserine", "linezolid", "bedaquiline", "clofazimine", "delamanid"]  

table_a = pandas.read_csv(snakemake.input[0] , index_col=0, sep="\t", header=0)

for sample, row in table_a.iterrows():
    for drug in drug_list:
        if table_a.at[sample, drug] == '-':
            table_a.at[sample, drug] = 'S'
        else:
            table_a.at[sample, drug] = 'R'

table_a.to_csv(snakemake.output[0], sep="\t")
