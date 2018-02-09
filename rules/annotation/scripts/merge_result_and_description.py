import pandas

res = pandas.read_csv(snakemake.input[0], delimiter="\t")
db = pandas.read_csv(snakemake.input[1], delimiter="\t")

db2 = pandas.DataFrame({"gene_uniprot":db[["gene", "uniprot_accession"]].apply(lambda x: '_'.join(x), axis=1), "description":db["description"], "gene":db["gene"]})

merged = pandas.merge(res, db2, left_on="virulence_factor_ID", right_on="gene_uniprot")


cols = list(merged)

cols.insert(1, cols.pop(cols.index("description")))

merged = merged.ix[:, cols]
merged.to_csv(snakemake.output[0], header=True, sep="\t", index=False)

writer = pandas.ExcelWriter(snakemake.output[1])
merged.to_excel(writer, snakemake.wildcards["sample"], index=False)
writer.save()
