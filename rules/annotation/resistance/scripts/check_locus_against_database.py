from Bio import SeqIO, Seq
import pandas

annotations = pandas.read_csv(snakemake.input["db"], sep="\t")
#the DB of annotation, see examples for file specifications. Compulsory columns are Gene, Postion, WildTypeAminoAcidOrNucleotide, MutationType
#Only works for single chromosome reference

gene_to_locustag = pandas.read_csv(snakemake.input["gene_to_locus"], sep="\t", index_col=0)

gbk_features = {}

incorrect_column_names = ["Gene", "Position", "WildTypeAminoAcidOrNucleotide", "MutatedAminoAcidOrNucleotide", "MutationType"]

incorrect_annotations = pandas.DataFrame(columns=incorrect_column_names)

bed_column_names = ["Reference", "Start", "End", "Gene", "Strand", "PositionInGene", "ReferenceAminoAcidOrNucleotide", "ReferenceCodonOrNucleotide"]

bed_file = pandas.DataFrame(columns=bed_column_names)

#We create a dictionnary with locus tag as key to access our genome of reference
#This part should be changed for references with many chromosmes

for record in SeqIO.parse(snakemake.input["gbk"], "genbank"):
    genome_reference = record
    version = genome_reference.id
    for feat in record.features:
        if feat.type=="gene":
            gbk_features[feat.qualifiers["locus_tag"][0]] = feat

#For each annotation, we check that the position in the genome corresponds to the annotated wildtype. We can thus verify that we know the coordinates in the reference for this position. There are six different cases, a mutation in the coding sequence of a gene (codon), in the promoter (nucleotide) or in a ribosomal gene (nucleotide), which can be on the positive or negative strand. Numbering scheme should be: 1 for first nucleotide of the start codon, and the previous nucleotide is -1. Otherwise this checking will fail. If the positions match, we construct the bed file with respect to the reference genome. 

for gene in annotations["Gene"].unique():
    for index, row in annotations.loc[(annotations["Gene"] == gene) & (annotations["MutationType"] == "SNP")].iterrows():
        status = 0
        #Verify that locus tag is present in genome
        if gene_to_locustag.loc[gene, "LocusTag"] in gbk_features.keys():
            coordinates = gbk_features[gene_to_locustag.loc[gene, "LocusTag"]].location
            position_in_gene = int(row["Position"])
            #RNA gene checking: nucleotide and not in promoter (although very unlikely to be annotated)
            if "rrs" in gene and position_in_gene > 0:
                wildtype = row["WildTypeAminoAcidOrNucleotide"].lower()
                codon_or_nucleotide = wildtype
                if coordinates.strand == 1 and record.seq[coordinates.start + position_in_gene - 1].lower() == row["WildTypeAminoAcidOrNucleotide"].lower():
                    status = 1
                    coordinates_for_bed = (coordinates.start + position_in_gene -1, coordinates.start + position_in_gene, coordinates.strand)
                #This is tested with fake E coli data in the repo
                elif coordinates.strand == -1 and Seq.Seq(record.seq[coordinates.end - position_in_gene]).reverse_complement().lower() ==  row["WildTypeAminoAcidOrNucleotide"].lower():
                    status = 1
                    coordinates_for_bed = (coordinates.end - position_in_gene, coordinates.end - position_in_gene + 1, coordinates.strand)
            #Promoter regions checking: nucleotide, negative positions
            elif position_in_gene < 0:
                wildtype = row["WildTypeAminoAcidOrNucleotide"].lower()
                codon_or_nucleotide = wildtype
                if coordinates.strand == 1 and record.seq[coordinates.start + position_in_gene].lower() == row["WildTypeAminoAcidOrNucleotide"].lower():
                    status = 1
                    coordinates_for_bed = (coordinates.start + position_in_gene, coordinates.start + position_in_gene + 1, coordinates.strand)
                elif coordinates.strand == -1 and Seq.Seq(record.seq[coordinates.end - position_in_gene -1]).reverse_complement().lower() ==  row["WildTypeAminoAcidOrNucleotide"].lower():
                    status = 1
                    coordinates_for_bed = (coordinates.end - position_in_gene - 1, coordinates.end - position_in_gene, coordinates.strand)
            #The rest is codon
            else:
                wildtype = row["WildTypeAminoAcidOrNucleotide"].upper()
                if coordinates.strand == 1 and record.seq[(coordinates.start + (position_in_gene - 1)*3):(coordinates.start + position_in_gene*3)].translate(table=11).lower() == row["WildTypeAminoAcidOrNucleotide"].lower():
                    status = 1
                    coordinates_for_bed = (coordinates.start + (position_in_gene - 1)*3, coordinates.start + position_in_gene*3, coordinates.strand)
                    codon_or_nucleotide = str(record.seq[(coordinates.start + (position_in_gene - 1)*3):(coordinates.start + position_in_gene*3)])
                elif coordinates.strand == -1 and record.seq[(coordinates.end - position_in_gene * 3):(coordinates.end - (position_in_gene - 1)*3)].reverse_complement().translate(table=11).lower() == row["WildTypeAminoAcidOrNucleotide"].lower():
                    status = 1
                    coordinates_for_bed = (coordinates.end - position_in_gene*3, coordinates.end - (position_in_gene - 1)*3, coordinates.strand)
                    codon_or_nucleotide = str(record.seq[(coordinates.end - position_in_gene * 3):(coordinates.end - (position_in_gene - 1)*3)].reverse_complement())
        if status:
            bed_file = bed_file.append(pandas.DataFrame([[version, coordinates_for_bed[0], coordinates_for_bed[1], gene, int(coordinates_for_bed[2]), row["Position"], wildtype, codon_or_nucleotide]], columns=bed_column_names))
        else:
            incorrect_annotations = incorrect_annotations.append(pandas.DataFrame([[gene, row["Position"], row["WildTypeAminoAcidOrNucleotide"], row["MutatedAminoAcidOrNucleotide"], row["MutationType"]]], columns=incorrect_column_names))

bed_file.sort_values(by=["Reference", "Start"], ascending=[True, True]).drop_duplicates().to_csv(snakemake.output["bed_correct"], sep="\t", index=False, header=False)
incorrect_annotations.to_csv(snakemake.output["incorrect_annotation"], sep="\t", index=False)
