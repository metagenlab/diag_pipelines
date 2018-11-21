
import numpy

fna_file = snakemake.input[0]
gbk_file = snakemake.input[1]
depth_file = snakemake.input[2]

out_CDS_depth = snakemake.output[0]
out_contig_depth = snakemake.output[1]

def get_contig_name2contig_coverage(samtool_depth_file):
    contig2depth = {}
    for line in open(samtool_depth_file, 'r'):
        data = line.rstrip().split('\t')
        if data[0] not in contig2depth:
            contig2depth[data[0]] = [int(data[2])]
        else:
            contig2depth[data[0]].append(int(data[2]))

    return contig2depth

def get_contig_id2median_depth(contig2depth_list):
    contig2med = {}
    for contig in contig2depth_list:
        m = numpy.median(contig2depth_list[contig])
        contig2med[contig] = m
    compl_list = numpy.concatenate(list(contig2depth_list.values()))
    contig2med["TOTAL"] = round(numpy.median(compl_list), 2)
    return contig2med

def get_contig_id2mean_depth(contig2depth_list):
    contig2mean = {}
    for contig in contig2depth_list:
        m = numpy.mean(contig2depth_list[contig])
        contig2mean[contig] = m
    compl_list = numpy.concatenate(list(contig2depth_list.values()))

    contig2mean["TOTAL"] = round(numpy.mean(compl_list), 2)
    return contig2mean

def contig_id2contig_length(contigs_file):

    from Bio import SeqIO

    fasta_handle = open(contigs_file, 'r')
    id2length = {}

    for record in SeqIO.parse(fasta_handle, "fasta"):
        id2length[record.name] = len(record)
    id2length["TOTAL"] = sum(list(id2length.values()))
    return id2length

def contig_id2gc_content(contigs_file):

    from Bio import SeqIO
    from Bio.SeqUtils import GC

    fasta_handle = open(contigs_file, 'r')
    id2gc = {}
    cumul_seq = ''
    for record in SeqIO.parse(fasta_handle, "fasta"):
        id2gc[record.name] = GC(record.seq)
        cumul_seq += record.seq
    id2gc["TOTAL"] = round(GC(cumul_seq), 2)
    return id2gc

def gbk2CDS_loc(gbk_file):

    from Bio import SeqIO
    record2locus2location = {}
    with open(gbk_file, 'r') as f:
        for record in SeqIO.parse(f, 'genbank'):
            record2locus2location[record.name] = {}
            for feature in record.features:
                if feature.type == 'CDS':
                    startstop = [feature.location.start, feature.location.end]
                    locus_tag = feature.qualifiers['locus_tag'][0]
                    record2locus2location[record.name][locus_tag] = startstop
    return record2locus2location

def parse_smatools_depth(samtools_depth):
    import pandas

    with open(samtools_depth, 'r') as f:
        table = pandas.read_csv(f, sep='\t', header=None, index_col=0)
    return table


### MAIN ###
contig_id2contig_length = contig_id2contig_length(fna_file)
contig_id2contig_gc = contig_id2gc_content(fna_file)
contig_id2depth = get_contig_name2contig_coverage(depth_file)
contig2median_depth = get_contig_id2median_depth(contig_id2depth)
contig2mean_depth = get_contig_id2mean_depth(contig_id2depth)
record2gene2coord = gbk2CDS_loc(gbk_file)

samtools_dataframe = parse_smatools_depth(depth_file)
median_depth = numpy.median(samtools_dataframe.iloc[:, 1])

# write gene depth
with open(out_CDS_depth, 'w') as f:
    f.write("contig\tgene\tstart\tend\tdepth\tratio_assembly\tcontig_depth\tcontig_ratio_depth\tcontig_length\n")
    for record in record2gene2coord:
        # get record subset
        subset_table = samtools_dataframe.loc[record]

        for gene in record2gene2coord[record]:
            # attention range
            # index starts from 0
            # range does not include "last count"
            start_pos = int(record2gene2coord[record][gene][0]) - 1
            end_pos = int(record2gene2coord[record][gene][1])
            if start_pos > end_pos:
                start_pos = int(record2gene2coord[record][gene][1]) - 1
                end_pos = int(record2gene2coord[record][gene][0])

            gene_median = numpy.median(subset_table.iloc[start_pos:end_pos, 1])
            f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
            record, gene, start_pos, end_pos, gene_median, round(gene_median / median_depth, 2), contig2median_depth[record],
            round(contig2median_depth[record] / median_depth, 2), contig_id2contig_length[record]))

# write contig depth
with open(out_contig_depth, 'w') as g:
    g.write("contig\tmean_depth\tmedian_depth\tgc_content\tcontig_size\n")
    # only write data for filtered contigs
    for contig in record2gene2coord:
        g.write("%s\t%s\t%s\t%s\t%s\n" % (contig,
                                  contig2mean_depth[contig],
                                  contig2median_depth[contig],
                                  contig_id2contig_gc[contig],
                                  contig_id2contig_length[contig]))

    g.write("%s\t%s\t%s\t%s\t%s\n" % ("TOTAL",
                              contig2mean_depth["TOTAL"],
                              contig2median_depth["TOTAL"],
                              contig_id2contig_gc["TOTAL"],
                              contig_id2contig_length["TOTAL"]))
