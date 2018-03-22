from Bio.Seq import Seq
import pandas

locus_tags = pandas.read_csv(snakemake.input['locus_tag'], sep="\t", index_col=1)


codons = pandas.read_csv(snakemake.input["resistance_codons"], sep="\t", header=None, names=["Annotation", "Reference", "Alternative"])

not_codons = pandas.read_csv(snakemake.input["resistance_not_codons"], sep="\t", header=None, names=["Annotation", "NuclReference", "NuclAlternative"])


codons = codons.assign(LocusTag= [str(x.split(":")[0].replace(">", "")) for x in codons["Annotation"]])
codons = codons.assign(CodonPosition= [str(x.split(":")[1]) for x in codons["Annotation"]])

not_codons = not_codons.assign(LocusTag= [str(x.split(":")[0].replace(">", "")) for x in not_codons["Annotation"]])
not_codons = not_codons.assign(NuclPosition= [int(x.split(":")[1]) for x in not_codons["Annotation"]])


codons = codons.assign(ReferenceAA= [str(Seq(x).translate()) for x in codons["Reference"]])
codons = codons.assign(AlternativeAA= [str(Seq(x).translate()) for x in codons["Alternative"]])
codons = codons.assign(AAPosition = [int(int(x.split("-")[1])//3) for x in codons["CodonPosition"]])

codons = codons.set_index("LocusTag")
not_codons = not_codons.set_index("LocusTag")

writer = pandas.ExcelWriter(snakemake.output["resistance"])

pandas.concat([codons, not_codons]).join(locus_tags)[["Gene", "CodonPosition", "ReferenceAA", "AAPosition", "AlternativeAA", "NuclReference", "NuclPosition", "NuclAlternative"]].to_excel(writer, index=True)

writer.save()
