
import pandas 

df_list = [pandas.read_csv(table, sep="\t", index_col=None, header=0) for table in snakemake.input]

combined = result = pandas.concat(df_list)

combined_filtered = combined.drop_duplicates('ORF_ID', keep='first')

combined_filtered.to_csv(snakemake.output[0], sep="\t", index=None)