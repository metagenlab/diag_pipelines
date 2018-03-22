from Bio.Seq import Seq
import pandas

locus_tags = pandas.read_csv(snakemake.input['locus_tag'], sep="\t", index_col=1)

nucleotides = pandas.read_excel(snakemake.input["formated_nucleotides"], sheet_name=snakemake.wildcards["sample"], header=0, index_col=0)
codons = pandas.read_excel(snakemake.input["formated_aa"], sheet_name=snakemake.wildcards["sample"], header=0, index_col=0)

writer = pandas.ExcelWriter(snakemake.output["resistance"])

pandas.concat([nucleotides, codons]).join(locus_tags)[["Gene", "CodonPosition", "ReferenceAA", "AAPosition", "AlternativeAA", "NuclReference", "NuclPosition", "NuclAlternative"]].to_excel(writer, index=True, sheet_name=snakemake.wildcards["sample"])

writer.save()
