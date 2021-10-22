
import pandas
from Bio import SeqIO

if "flash_hist" in snakemake.input:   
    flash_hist = snakemake.input["flash_hist"] #samples/{sample}/reads/raw/{sample}.
else:
    flash_hist = []
contigs_files = snakemake.input["contigs"] # samples/{sample}/assembly/spades/contigs.fasta
contigs_depth = snakemake.input["contigs_depth"] 
mash_results = snakemake.input["mash_results"] # samples/{sample}/contamination/mash/assembly/distances_formated_no_virus.tsv
centrifuge_tables = snakemake.input["centrifuge_tables"] # report/contamination/centrifuge/{sample}/centrifuge_kraken_format.txt

o = open(snakemake.output[0], "w")

sample2max_flash = {}
sample2centrifuge = {}
sample2mash = {}
sample2small_contigs = {}
contig2contigs_depth = {}

# parse flash histograms to extract max of the distribution
for flash in flash_hist:
    sample = flash.split("/")[1]
    # structure: n, fragment_length, abundance
    table = pandas.read_csv(flash, sep="\t", header=None)
    # extract the most abundant length
    sample2max_flash[sample] = table[0][table[1].idxmax()]

# parse contig depth 
#contig_1	759142	60952585	80.29141451796897	22.849062292288043
for contig in contigs_depth:
    sample = contig.split("/")[1]
    contig2contigs_depth[sample] = {}
    contig2contigs_depth[sample]["lower_5"] = 0
    contig2contigs_depth[sample]["lower_10"] = 0
    contig2contigs_depth[sample]["lower_15"] = 0

    table = pandas.read_csv(contig, sep="\t", header=None)
    depth_largest_contig = table.iloc[0,3]
    contig2contigs_depth[sample]["depth_largest"] = depth_largest_contig

    for n,row in table.iterrows():
        contig_name = row[0]
        depth = float(row[3])
        if depth < 5:
            contig2contigs_depth[sample]["lower_5"] += 1
        if depth < 10:
            contig2contigs_depth[sample]["lower_10"] += 1
        if depth < 15:
            contig2contigs_depth[sample]["lower_15"] += 1

# parse contig length
for contigs in contigs_files:
    sample = contigs.split("/")[1]
    sample2small_contigs[sample] = {}
    sample2small_contigs[sample]["smaller_500"] = 0
    sample2small_contigs[sample]["smaller_1000"] = 0
    sample2small_contigs[sample]["smaller_1500"] = 0
    sample2small_contigs[sample]["smaller_2000"] = 0
    sample2small_contigs[sample]["total"] = 0

    records = SeqIO.parse(contigs, "fasta")

    for record in records:
        sample2small_contigs[sample]["total"] += 1
        l = len(record.seq)
        if l < 500:
            sample2small_contigs[sample]["smaller_500"] += 1
        if l < 1000:
            sample2small_contigs[sample]["smaller_1000"] += 1
        if l < 1500:
            sample2small_contigs[sample]["smaller_1500"] += 1
        if l < 2000:
            sample2small_contigs[sample]["smaller_2000"] += 1


# parse mash results
for mash in mash_results:
    sample = mash.split("/")[1]
    sample2mash[sample] = {}
    # 0.999281	985/1000	0	Pseudomonas aeruginosa strain AZPAE14893 AZPAE14893_contig_1, whole genome shotgun sequence
    # identity, shared-hashes,
    table = pandas.read_csv(mash, sep="\t", header=None)
    best_hit_identity = table.iloc[0,0]
    best_hit_shared_hash = table.iloc[0,1]
    best_hit_description = " ".join(table.iloc[0,3].split(" ")[0:2])
    sample2mash[sample]["best_hit_identity"] = best_hit_identity
    sample2mash[sample]["best_hit_shared_hash"] = best_hit_shared_hash.split("/")[0]
    sample2mash[sample]["best_hit_description"] = best_hit_description

# parse centrifuge results
for centrifuge in centrifuge_tables:
    sample = centrifuge.split("/")[1]
    sample2centrifuge[sample] = {}
    table = pandas.read_csv(centrifuge, sep="\t")
    # calculate ratio of unique reads for each row
    table["proportion"] = (table["numUniqueReads"]/sum(table["numUniqueReads"]))*100
    # sort by proportion to extract the 3 hits with highest abundance
    table = table.sort_values(by=['proportion'], ascending=False)
    # species name: .iloc[0,0]
    # abunance: .iloc[0,7]
    sample2centrifuge[sample]["hit_1"] = [table.iloc[0,0], table.iloc[0,7]]
    sample2centrifuge[sample]["hit_2"] = [table.iloc[1,0], table.iloc[1,7]]
    sample2centrifuge[sample]["hit_3"] = [table.iloc[2,0], table.iloc[2,7]]



header = [
          "sample",
          "centrifuge_hit_1_name",
          "centrifuge_hit_1_proportion",
          "centrifuge_hit_2_name",
          "centrifuge_hit_2_proportion",
          "centrifuge_hit_3_name",
          "centrifuge_hit_3_proportion",
          "mash_best_hit_identity",
          "mash_best_hit_shared_hash",
          "mash_best_hit_description",
          "n_contigs",
          "contigs_smaller_500",
          "contigs_smaller_1000",
          "contigs_smaller_2000",
          "depth_lower_5",
          "depth_lower_10",
          "depth_lower_15"
          ]

if "flash_hist" in snakemake.input:  
    header.append("max_flash")

o.write('\t'.join(header) + "\n")
for sample in sample2centrifuge:

    # centrifuge
    centrifuge_hit_1_name, centrifuge_hit_1_proportion = sample2centrifuge[sample]["hit_1"] 
    centrifuge_hit_2_name, centrifuge_hit_2_proportion = sample2centrifuge[sample]["hit_2"] 
    centrifuge_hit_3_name, centrifuge_hit_3_proportion = sample2centrifuge[sample]["hit_3"] 
    # mash
    mash_best_hit_identity = sample2mash[sample]["best_hit_identity"]
    mash_best_hit_shared_hash = sample2mash[sample]["best_hit_shared_hash"]
    mash_best_hit_description = sample2mash[sample]["best_hit_description"]
    # small contigs
    n_contigs = sample2small_contigs[sample]["total"]
    contigs_smaller_500 = sample2small_contigs[sample]["smaller_500"]
    contigs_smaller_1000 = sample2small_contigs[sample]["smaller_1000"]
    contigs_smaller_1500 = sample2small_contigs[sample]["smaller_1500"]
    contigs_smaller_2000 = sample2small_contigs[sample]["smaller_2000"]
    # depth
    depth_lower_5 = contig2contigs_depth[sample]["lower_5"]
    depth_lower_10 = contig2contigs_depth[sample]["lower_10"]
    depth_lower_15 = contig2contigs_depth[sample]["lower_15"]

    row = [
          sample,
          centrifuge_hit_1_name,
          centrifuge_hit_1_proportion,
          centrifuge_hit_2_name,
          centrifuge_hit_2_proportion,
          centrifuge_hit_3_name,
          centrifuge_hit_3_proportion,
          mash_best_hit_identity,
          mash_best_hit_shared_hash,
          mash_best_hit_description,
          n_contigs,
          contigs_smaller_500,
          contigs_smaller_1000,
          contigs_smaller_2000,
          depth_lower_5,
          depth_lower_10,
          depth_lower_15
    ]
    if "flash_hist" in snakemake.input:
        # flash
        max_flash = sample2max_flash[sample]
        row.append(max_flash) 
    o.write('\t'.join([str(i) for i in row]) + "\n")
o.close() 
