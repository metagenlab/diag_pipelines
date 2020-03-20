
import pandas 

COLUMNS = ["isoniazid","rifampicin","ethambutol","pyrazinamide"]

table_ref = pandas.read_csv(snakemake.params["reference_table"] , index_col=0, sep="\t")

o = open(snakemake.output[0], "w")

with open(snakemake.input[0], "r") as f:
    for row in f:
        tool, sample, drug, gene, position, change, vartype = row.rstrip().split("\t")
        if drug in COLUMNS:
            phenotype = table_ref.at[sample, drug]
            o.write(f"{tool}\t{sample}\t{drug}\t{gene}\t{position}\t{change}\t{vartype}\t{phenotype}\n")
        else:
            o.write(f"{tool}\t{sample}\t{drug}\t{gene}\t{position}\t{change}\t{vartype}\t-\n")
