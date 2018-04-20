from Bio import SeqIO
import pandas

all_features = {}
for record in SeqIO.parse(snakemake.input["gbk"], "genbank"):
    h37rv = record
    for feat in record.features:
        if feat.type=="gene":
            all_features[feat.qualifiers["locus_tag"][0]] = feat


mutations = pandas.read_csv(snakemake.input["db_correct"], sep="\t", index_col=0)
correspondance = pandas.read_csv(snakemake.input["gene_to_locus"], sep="\t", index_col=0)

print(mutations)



for index, row in mutations.join(correspondance).iterrows():
    coordinates = all_features[row["LocusTag"]].location
    if index == "rrs":
        if coordinates.strand == 1:
            loc = (coordinates.start + row["PositionMTB"], coordinates.start + row["PositionMTB"] + 1)
        elif coordinates.strand == -1:
            loc = (coordinates.end - row["PositionMTB"], coordinates.end - row["PositionMTB"] - 1)
        print(loc)
    elif int(row["PositionMTB"]) < 0:
        if coordinates.strand == 1:
            loc = (coordinates.start + row["PositionMTB"], coordinates.start + row["PositionMTB"] + 1)
        elif coordinates.strand == -1:
            loc = (coordinates.end - row["PositionMTB"], coordinates.end - row["PositionMTB"] + 1)
        print(loc)
    else:
        if coordinates.strand == 1:
            loc = (coordinates.start + (row["PositionMTB"]-1)*3, coordinates.start + row["PositionMTB"]*3)
        elif coordinates.strand == -1:
            loc = (coordinates.end - row["PositionMTB"]*3, coordinates.end - (row["PositionMTB"]-1)*3)
    print(row, loc)
    


