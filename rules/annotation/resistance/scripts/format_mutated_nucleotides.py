import pandas

not_codons = pandas.read_csv(snakemake.input["resistance_nucleotides"], sep="\t", header=None, names=["Annotation", "NuclReference", "NuclAlternative"])
not_codons = not_codons.assign(LocusTag= [str(x.split(":")[0].replace(">", "")) for x in not_codons["Annotation"]])
not_codons = not_codons.assign(NuclPosition= [int(x.split(":")[1].split("-")[1]) for x in not_codons["Annotation"]])
not_codons = not_codons.set_index("LocusTag")

writer = pandas.ExcelWriter(snakemake.output["formated_nucleotides"])

not_codons[["NuclReference", "NuclPosition", "NuclAlternative"]].to_excel(writer, index=True, sheet_name=snakemake.wildcards["sample"])

writer.save()
