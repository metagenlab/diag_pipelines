import pandas


writer = pandas.ExcelWriter(snakemake.output["xlsx"])
pandas.read_csv(snakemake.input["tsv"], sep="\t").to_excel(writer, index=False, sheet_name=snakemake.wildcards["software"])
writer.save()
