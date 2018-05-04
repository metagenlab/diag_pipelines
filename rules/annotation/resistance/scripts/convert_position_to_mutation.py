from Bio import SeqIO, Seq

import pandas

#We read the genotype info, which is always per nucleotide and not per codon. If we have more than one mutation on the same codon, we annotate all mutated positions. Genotype however can be codon or nucleotide.

genotype = pandas.read_csv(snakemake.input["genotype"], sep="\t", names=["Chrom", "Position", "Genotype", "Strand"])

aa_or_nucl = []
codon_or_nucl = []

#If the length of the Genotype is 3, we search for the correct orientation of the codon and translate it to get the mutated AA. We also annotate with the correct value of the codon in the context of the gene, which will be different from the VCF ALT values for -1 strand genes.

for index, row in genotype.iterrows():
    if len(row["Genotype"])==3:
        if row["Strand"]==1:
            aa_or_nucl.append(str(Seq.Seq(row["Genotype"]).translate(table=11)))
            codon_or_nucl.append(str(Seq.Seq(row["Genotype"])))
        elif row["Strand"]==-1:
            aa_or_nucl.append(str(Seq.Seq(row["Genotype"]).reverse_complement().translate(table=11)))
            codon_or_nucl.append(str(Seq.Seq(row["Genotype"]).reverse_complement()))
    elif len(row["Genotype"])==1:
        if row["Strand"]==1:
            aa_or_nucl.append(row["Genotype"])
            codon_or_nucl.append(row["Genotype"])
        elif row["Strand"]==-1:
            aa_or_nucl.append(str(Seq.Seq(row["Genotype"]).reverse_complement()))
            codon_or_nucl.append(str(Seq.Seq(row["Genotype"]).reverse_complement()))

#We switch to bed file coordinates, which are 0-based, whereas previous vcf coordinates are 1-based.

genotype["Start"] = genotype["Position"] - 1
genotype["End"] = genotype["Position"]

pandas.concat([genotype[["Chrom", "Start", "End"]], pandas.Series(aa_or_nucl), pandas.Series(codon_or_nucl)], axis=1).to_csv(snakemake.output["mutation"], sep="\t", index=False, header=False)

