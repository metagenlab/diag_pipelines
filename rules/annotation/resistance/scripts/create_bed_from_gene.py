from Bio import SeqIO
from Bio import Seq
import pandas

gene_to_locustag = pandas.read_csv(snakemake.input["gene_to_locus"], sep="\t", index_col=0)
gbk_features = {}

for record in SeqIO.parse(snakemake.input["gbk"], "genbank"):
    h37rv = record
    version = h37rv.id
    for feat in record.features:
        if feat.type=="gene":
            gbk_features[feat.qualifiers["locus_tag"][0]] = feat

coordinates = gbk_features[gene_to_locustag.loc[snakemake.wildcards["gene"], "LocusTag"]].location

with open(snakemake.output["bed"], "w") as f:
    f.write(version + "\t" + str(coordinates.start - snakemake.params["up_down"]) + "\t" + str(coordinates.end + snakemake.params["up_down"]) + "\t" + str(coordinates.strand) + "\n")
