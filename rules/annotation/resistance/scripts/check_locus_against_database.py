from Bio import SeqIO
from Bio import Seq
import pandas

mutation_entries = pandas.read_csv(snakemake.input["db"], sep="\t")
correspondance = pandas.read_csv(snakemake.input["gene_to_locus"], sep="\t", index_col=0)

all_features = {}

column_names = ["Gene", "PositionMTB", "WildTypeAminoAcidOrNucleotide", "MutatedAminoAcidOrNucleotide", "MutationType"]
correct = pandas.DataFrame(columns=column_names)
incorrect = pandas.DataFrame(columns=column_names)

for record in SeqIO.parse(snakemake.input["gbk"], "genbank"):
    h37rv = record
    for feat in record.features:
        if feat.type=="gene":
            all_features[feat.qualifiers["locus_tag"][0]] = feat

for gene in mutation_entries["Gene"].unique():
    if correspondance.loc[gene, "LocusTag"] not in all_features.keys():
        incorrect = incorrect.append(mutation_entries.loc[(mutation_entries["Gene"] == gene), ["Gene", "PositionMTB", "WildTypeAminoAcidOrNucleotide", "MutatedAminoAcidOrNucleotide", "MutationType"]])
    else:
        for index, row in mutation_entries.loc[(mutation_entries["Gene"] == gene) & (mutation_entries["MutationType"] == "SNP")].iterrows():
            status = 0
            coordinates = all_features[correspondance.loc[gene, "LocusTag"]].location
            #RNA gene checking
            if gene == "rrs":
                if coordinates.strand == 1 and record.seq[coordinates.start:coordinates.end][int(row["PositionMTB"])-1].lower() == row["WildTypeAminoAcidOrNucleotide"].lower():
                    status = 1
                elif coordinates.strand == -1 and Seq.Seq(record.seq[coordinates.end-int(row["PositionMTB"])-1]).reverse_complement().lower() ==  row["WildTypeAminoAcidOrNucleotide"].lower():
                    status = 1
            #Promoter regions checking
            elif int(row["PositionMTB"]) < 0:
                if coordinates.strand == 1 and record.seq[coordinates.start+int(row["PositionMTB"])].lower() == row["WildTypeAminoAcidOrNucleotide"].lower():
                    status = 1
                elif coordinates.strand == -1 and Seq.Seq(record.seq[coordinates.end-int(row["PositionMTB"])-1]).reverse_complement().lower() ==  row["WildTypeAminoAcidOrNucleotide"].lower():
                    status = 1
            #Mutation leading to codon change checking
            else:
                if coordinates.strand == 1 and record.seq[coordinates.start:coordinates.end].translate(table=11, cds=True)[int(row["PositionMTB"])-1].lower() == row["WildTypeAminoAcidOrNucleotide"].lower():
                    status = 1 
                elif coordinates.strand == -1 and record.seq[coordinates.start:coordinates.end].reverse_complement().translate(table=11, cds=True)[int(row["PositionMTB"])-1].lower() == row["WildTypeAminoAcidOrNucleotide"].lower():
                    status = 1
            if status:
                correct = correct.append(pandas.DataFrame([[gene, row["PositionMTB"], row["WildTypeAminoAcidOrNucleotide"], row["MutatedAminoAcidOrNucleotide"], row["MutationType"]]], columns=column_names))
            else:
                incorrect = incorrect.append(pandas.DataFrame([[gene, row["PositionMTB"], row["WildTypeAminoAcidOrNucleotide"], row["MutatedAminoAcidOrNucleotide"], row["MutationType"]]], columns=column_names))
                
        
                                    
correct.to_csv(snakemake.output["correct_annotation"], sep="\t", index=False)
incorrect.to_csv(snakemake.output["incorrect_annotation"], sep="\t", index=False)

                                

                    

        
