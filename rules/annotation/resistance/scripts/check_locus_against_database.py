from Bio import SeqIO
from Bio import Seq
import pandas

annotations = pandas.read_csv(snakemake.input["db"], sep="\t")
gene_to_locustag = pandas.read_csv(snakemake.input["gene_to_locus"], sep="\t", index_col=0)

gbk_features = {}

column_names = ["Gene", "Position", "WildTypeAminoAcidOrNucleotide", "MutatedAminoAcidOrNucleotide", "MutationType"]

incorrect_annotations = pandas.DataFrame(columns=column_names)
bed_file = pandas.DataFrame(columns=["Reference", "Start", "End", "Gene", "Strand", "PositionInGene", "ReferenceAminoAcidorNucleotide"])

for record in SeqIO.parse(snakemake.input["gbk"], "genbank"):
    h37rv = record
    version = h37rv.id
    for feat in record.features:
        if feat.type=="gene":
            gbk_features[feat.qualifiers["locus_tag"][0]] = feat

for gene in annotations["Gene"].unique():
    if gene_to_locustag.loc[gene, "LocusTag"] not in gbk_features.keys():
        incorrect_annotations = incorrect_annotations.append(annotations.loc[(mutation_entries["Gene"] == gene), ["Gene", "Position", "WildTypeAminoAcidOrNucleotide", "MutatedAminoAcidOrNucleotide", "MutationType"]])
    else:
        for index, row in annotations.loc[(annotations["Gene"] == gene) & (annotations["MutationType"] == "SNP")].iterrows():
            status = 0
            coordinates = gbk_features[gene_to_locustag.loc[gene, "LocusTag"]].location
            position_in_gene = int(row["Position"])
            #RNA gene checking: no codon
            if "rrs" in gene:
                wildtype = row["WildTypeAminoAcidOrNucleotide"].lower()
                if coordinates.strand == 1 and record.seq[coordinates.start + position_in_gene - 1].lower() == row["WildTypeAminoAcidOrNucleotide"].lower():
                    status = 1
                    coordinates_for_bed = (coordinates.start + position_in_gene -1, coordinates.start + position_in_gene, coordinates.strand)
                #THIS IS TESTED WITH E_COLI TEST DB IN THE REPO
                elif coordinates.strand == -1 and Seq.Seq(record.seq[coordinates.end - position_in_gene]).reverse_complement().lower() ==  row["WildTypeAminoAcidOrNucleotide"].lower():
                    status = 1
                    coordinates_for_bed = (coordinates.end - position_in_gene, coordinates.end - position_in_gene + 1, coordinates.strand)
            #Promoter regions checking: no codon, negative positions
            elif position_in_gene < 0:
                wildtype = row["WildTypeAminoAcidOrNucleotide"].lower()                  
                if coordinates.strand == 1 and record.seq[coordinates.start + position_in_gene ].lower() == row["WildTypeAminoAcidOrNucleotide"].lower():
                    status = 1
                    coordinates_for_bed = (coordinates.start + position_in_gene, coordinates.start + position_in_gene + 1, coordinates.strand)
                elif coordinates.strand == -1 and Seq.Seq(record.seq[coordinates.end - position_in_gene -1]).reverse_complement().lower() ==  row["WildTypeAminoAcidOrNucleotide"].lower():
                    status = 1
                    coordinates_for_bed = (coordinates.end - position_in_gene - 1, coordinates.end - position_in_gene, coordinates.strand)
            #The rest is codon checking
            else:
                wildtype = row["WildTypeAminoAcidOrNucleotide"].upper()
                if coordinates.strand == 1 and record.seq[(coordinates.start + (position_in_gene - 1)*3):(coordinates.start + position_in_gene*3)].translate(table=11).lower() == row["WildTypeAminoAcidOrNucleotide"].lower():
                    status = 1
                    coordinates_for_bed = (coordinates.start + (position_in_gene - 1)*3, coordinates.start + position_in_gene*3, coordinates.strand)
                elif coordinates.strand == -1 and record.seq[(coordinates.end - position_in_gene * 3):(coordinates.end - (position_in_gene - 1)*3)].reverse_complement().translate(table=11).lower() == row["WildTypeAminoAcidOrNucleotide"].lower():
                    status = 1
                    coordinates_for_bed = (coordinates.end - position_in_gene*3, coordinates.end - (position_in_gene - 1)*3, coordinates.strand)
            if status:
                bed_file = bed_file.append(pandas.DataFrame([[version, coordinates_for_bed[0], coordinates_for_bed[1], gene, int(coordinates_for_bed[2]), row["Position"], wildtype]], columns=["Reference", "Start", "End", "Gene", "Strand", "PositionInGene", "ReferenceAminoAcidorNucleotide"]))
            else:
                incorrect_annotations = incorrect_annotations.append(pandas.DataFrame([[gene, row["Position"], row["WildTypeAminoAcidOrNucleotide"], row["MutatedAminoAcidOrNucleotide"], row["MutationType"]]], columns=column_names))
                
        
                                    
bed_file.sort_values(by=["Start"]).to_csv(snakemake.output["bed_correct"], sep="\t", index=False, header=False)
incorrect_annotations.to_csv(snakemake.output["incorrect_annotation"], sep="\t", index=False)


                                

                    

        
