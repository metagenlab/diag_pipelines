from Bio import SeqIO, Seq

import pandas

genotype = pandas.read_csv(snakemake.input["genotype"], sep=" ", names=["Chrom", "Position", "Genotype", "Strand"])

aa_or_nucl = []
codon_or_nucl = []

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

#we switch to bed file coordinates, which are 0 based

genotype["Start"] = genotype["Position"] - 1
genotype["End"] = genotype["Position"]

pandas.concat([genotype[["Chrom", "Start", "End"]], pandas.Series(aa_or_nucl), pandas.Series(codon_or_nucl)], axis=1).to_csv(snakemake.output["mutation"], sep="\t", index=False, header=False)

