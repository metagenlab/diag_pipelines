

gbk_file = snakemake.input[0]
blast_file = snakemake.input[1]
VF_annotation = snakemake.input[2]
outfile = snakemake.output[0]
import pandas


def gbk2CDS_loc(gbk_file):

    from Bio import SeqIO
    locus2location = {}
    with open(gbk_file, 'r') as f:
        for record in SeqIO.parse(f, 'genbank'):
            locus2location[record.name] = {}
            for feature in record.features:
                if feature.type == 'CDS':
                    startstop = [record.name, feature.location.start, feature.location.end]
                    locus_tag = feature.qualifiers['locus_tag'][0]
                    locus2location[locus_tag] = startstop
    return locus2location

# VF accession: VFG048628
locus2VF_hit = pandas.read_csv(blast_file,
                           delimiter='\t',
                           header=0).set_index("ORF_ID").to_dict()["matching_sequence"]
# VFG037176
# VF0470
# plc
# Acinetobacter baumannii ACICU
# Phospholipase C
# phospholipase C
# YP_001844723

VF_accession2VF_id = pandas.read_csv(VF_annotation,
                           delimiter='\t',
                           names=["VFACCESSION","VFID","VF","SPECIES","DESCRIPTION","PRODUCT","SOURCE"]).set_index("VFACCESSION").to_dict()["VFID"]
VF_accession2VF_name = pandas.read_csv(VF_annotation,
                           delimiter='\t',
                           names=["VFACCESSION","VFID","VF","SPECIES","DESCRIPTION","PRODUCT","SOURCE"]).set_index("VFACCESSION").to_dict()["VF"]

locus2location = gbk2CDS_loc(gbk_file)

with open(outfile, 'w') as f:
    for locus in locus2VF_hit:
        vf_hit = locus2VF_hit[locus]
        label = "%s_%s" % (VF_accession2VF_id[vf_hit],
                            VF_accession2VF_name[vf_hit])
        f.write("%s\t%s\t%s\t%s\n" % (locus2location[locus][0],
                                   locus2location[locus][1],
                                   locus2location[locus][2],
                                   label))