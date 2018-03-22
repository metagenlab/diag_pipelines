from Bio.Seq import Seq
import pandas

codons = pandas.read_csv(snakemake.input["resistance_codons"], sep="\t", header=None, names=["Annotation", "Reference", "Alternative"])

codons = codons.assign(LocusTag= [str(x.split(":")[0].replace(">", "")) for x in codons["Annotation"]])
codons = codons.assign(CodonPosition= [str(x.split(":")[1]) for x in codons["Annotation"]])
codons = codons.assign(ReferenceAA= [str(Seq(x).translate()) for x in codons["Reference"]])
codons = codons.assign(AlternativeAA= [str(Seq(x).translate()) for x in codons["Alternative"]])
codons = codons.assign(AAPosition = [int(int(x.split("-")[1])//3) for x in codons["CodonPosition"]])
codons = codons.set_index("LocusTag")

writer = pandas.ExcelWriter(snakemake.output["formated_aa"])


codons.loc[codons["ReferenceAA"]!=codons["AlternativeAA"], ["CodonPosition", "ReferenceAA", "AAPosition", "AlternativeAA"]].to_excel(writer, index=True, sheet_name=snakemake.wildcards["sample"])

writer.save()
