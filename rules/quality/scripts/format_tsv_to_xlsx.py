import pandas


writer = pandas.ExcelWriter(snakemake.output["xlsx"])
pandas.read_csv(snakemake.input["tsv"], sep="\t", header=None, names=["Identity", "SharedHashes", "PValue", "Query"]).to_excel(writer, index=False, sheet_name=snakemake.wildcards["sample"])
writer.save()
