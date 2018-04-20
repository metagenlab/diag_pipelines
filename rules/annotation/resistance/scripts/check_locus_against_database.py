from Bio import SeqIO
from Bio import Seq
import pandas

annotations = pandas.read_csv(snakemake.input["db"], sep="\t")
gene_to_locustag = pandas.read_csv(snakemake.input["gene_to_locus"], sep="\t", index_col=0)

gbk_features = {}

column_names = ["Gene", "PositionMTB", "WildTypeAminoAcidOrNucleotide", "MutatedAminoAcidOrNucleotide", "MutationType"]

correct_annotations = pandas.DataFrame(columns=column_names)
incorrect_annotations = pandas.DataFrame(columns=column_names)

for record in SeqIO.parse(snakemake.input["gbk"], "genbank"):
    h37rv = record
    for feat in record.features:
        if feat.type=="gene":
            gbk_features[feat.qualifiers["locus_tag"][0]] = feat

for gene in annotations["Gene"].unique():
    if gene_to_locustag.loc[gene, "LocusTag"] not in gbk_features.keys():
        incorrect_annotations = incorrect_annotations.append(annotations.loc[(mutation_entries["Gene"] == gene), ["Gene", "PositionMTB", "WildTypeAminoAcidOrNucleotide", "MutatedAminoAcidOrNucleotide", "MutationType"]])
    else:
        for index, row in annotations.loc[(annotations["Gene"] == gene) & (annotations["MutationType"] == "SNP")].iterrows():
            status = 0
            coordinates = gbk_features[gene_to_locustag.loc[gene, "LocusTag"]].location
            #RNA gene checking: no codon
            position_in_gene = int(row["PositionMTB"])
            if "rrs" in gene:
                if coordinates.strand == 1 and record.seq[coordinates.start + position_in_gene - 1].lower() == row["WildTypeAminoAcidOrNucleotide"].lower():
                    status = 1
                #THIS IS TESTED WITH E_COLI TEST DB IN THE REPO
                elif coordinates.strand == -1 and Seq.Seq(record.seq[coordinates.end - position_in_gene]).reverse_complement().lower() ==  row["WildTypeAminoAcidOrNucleotide"].lower():
                    status = 1
            #Promoter regions checking: no codon, negative positions
            elif position_in_gene < 0:
                if coordinates.strand == 1 and record.seq[coordinates.start + position_in_gene ].lower() == row["WildTypeAminoAcidOrNucleotide"].lower():
                    status = 1
                elif coordinates.strand == -1 and Seq.Seq(record.seq[coordinates.end - position_in_gene -1]).reverse_complement().lower() ==  row["WildTypeAminoAcidOrNucleotide"].lower():
                    status = 1
            #The rest is codon checking
            else:
                if coordinates.strand == 1 and record.seq[(coordinates.start + (position_in_gene - 1)*3):(coordinates.start+(position_in_gene * 3))].translate(table=11).lower() == row["WildTypeAminoAcidOrNucleotide"].lower():
                    status = 1 
                elif coordinates.strand == -1 and record.seq[(coordinates.end - position_in_gene * 3):(coordinates.end - (position_in_gene - 1)*3)].reverse_complement().translate(table=11).lower() == row["WildTypeAminoAcidOrNucleotide"].lower():
                    status = 1
            if status:
                correct_annotations = correct_annotations.append(pandas.DataFrame([[gene, row["PositionMTB"], row["WildTypeAminoAcidOrNucleotide"], row["MutatedAminoAcidOrNucleotide"], row["MutationType"]]], columns=column_names))
            else:
                incorrect_annotations = incorrect_annotations.append(pandas.DataFrame([[gene, row["PositionMTB"], row["WildTypeAminoAcidOrNucleotide"], row["MutatedAminoAcidOrNucleotide"], row["MutationType"]]], columns=column_names))
                
        
                                    
correct_annotations.to_csv(snakemake.output["correct_annotation"], sep="\t", index=False)
incorrect_annotations.to_csv(snakemake.output["incorrect_annotation"], sep="\t", index=False)

                                

                    

        
