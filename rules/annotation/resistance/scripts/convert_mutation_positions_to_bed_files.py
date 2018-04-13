import pandas

mutations = pandas.read_csv(snakemake.input["db_correct"], sep="\t")
correspondance = pandas.read_csv(snakemake.input["gene_to_locus"], sep="\t", index_col=0)

codons = mutations.loc[(mutations["Gene"]!="rrs") & (mutations["PositionMTB"]>0)]
codons = codons.assign(Start= lambda x: (x.PositionMTB - 1)*3 + snakemake.params["upstream_downstream_size"])
codons = codons.assign(End= lambda x: (x.PositionMTB*3) + snakemake.params["upstream_downstream_size"])

rna = mutations.loc[(mutations["Gene"]=="rrs")]
rna = rna.assign(Start = lambda x: (x.PositionMTB - 1 + snakemake.params["upstream_downstream_size"]))
rna = rna.assign(End = lambda x: (x.PositionMTB + snakemake.params["upstream_downstream_size"]))

promoters = mutations.loc[(mutations["PositionMTB"]<0)]
promoters = promoters.assign(Start = lambda x: (x.PositionMTB + snakemake.params["upstream_downstream_size"]))
promoters = promoters.assign(End = lambda x: (x.PositionMTB + 1 + snakemake.params["upstream_downstream_size"]))

codons_nucleotides = codons.append(rna).append(promoters)
codons_nucleotides = codons_nucleotides.set_index("Gene")

codons_nucleotides = codons_nucleotides.join(correspondance)

codons_nucleotides[["LocusTag", "Start", "End"]].to_csv(snakemake.output["bed_shifted"], sep="\t", index=False, header=False)
