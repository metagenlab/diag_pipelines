import pandas


writer = pandas.ExcelWriter(snakemake.output["xlsx"])
pandas.read_csv(snakemake.input["tsv"], sep="\t").to_excel(writer, index=False, sheet_name=snakemake.wildcards["sample"]+"_"+snakemake.wildcards["software"])
writer.save()
