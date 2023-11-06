
import pandas 

df_list = []
for table in snakemake.input:
    df = pandas.read_csv(table, sep="\t", index_col=None, header=0)
    df["sample"] = table.split("/")[1]
    df_list.append(df)

combined = result = pandas.concat(df_list)

combined.to_csv(snakemake.output[0], sep="\t", index=None)