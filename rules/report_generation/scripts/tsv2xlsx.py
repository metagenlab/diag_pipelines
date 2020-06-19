
import pandas 


input = snakemake.input[0]
output = snakemake.output[0]

df = pandas.read_csv(input, sep="\t", header=0)
df.to_excel(output)
