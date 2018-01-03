import pandas

res = pandas.read_csv(snakemake.input[0], delimiter="\t")
db = pandas.read_csv(snakemake.input[1], delimiter="\t")

merged = pandas.merge(res, db, left_on="Virulence factor ID", right_on="gene")


merged[["Virulence factor ID", "Sequenced gene ID", "Amino acid percentage identity", "Coverage of the Virulence factor protein on the alignement", "E-value", "uniprot_accession", "description"]].to_csv(snakemake.output[0], header=True, sep="\t", index=False)

writer = pandas.ExcelWriter(snakemake.output[1])
merged[["Virulence factor ID", "Sequenced gene ID", "Amino acid percentage identity", "Coverage of the Virulence factor protein on the alignement", "E-value", "uniprot_accession", "description"]].to_excel(writer, snakemake.wildcards["sample"], index=False)
writer.save()
