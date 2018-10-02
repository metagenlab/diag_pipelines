
import numpy

fna_file = snakemake.input[0]
gbk_file = snakemake.input[1]
depth_file = snakemake.input[2]

out_CDS_depth = snakemake.output[0]
#out_contig_depth = snakemake.output[1]

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
    return contig2med

def contig_id2contig_length(contigs_file):

    from Bio import SeqIO

    fasta_handle = open(contigs_file, 'r')
    id2length = {}
    for record in SeqIO.parse(fasta_handle, "fasta"):
        id2length[record.name] = len(record)
    return id2length

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
contig_id2depth = get_contig_name2contig_coverage(depth_file)
contig2median_depth = get_contig_id2median_depth(contig_id2depth)
record2gene2coord = gbk2CDS_loc(gbk_file)

samtools_dataframe = parse_smatools_depth(depth_file)
median_depth = numpy.median(samtools_dataframe.iloc[:, 1])

with open(out_CDS_depth, 'w') as f:
    f.write("Contig\tgene\tstart\tend\tdepth\tratio_assembly\tcontig_depth\tcontig_ratio_depth\tcontig_length\n")
    for record in record2gene2coord:
        # get record subset
        subset_table = samtools_dataframe.loc[record]

        for gene in record2gene2coord[record]:
            # attention range
            # index commence a 0
            # range n'inclu pas le dernier chiffre
            start_pos = int(record2gene2coord[record][gene][0]) - 1
            end_pos = int(record2gene2coord[record][gene][1])
            if start_pos > end_pos:
                start_pos = int(record2gene2coord[record][gene][1]) - 1
                end_pos = int(record2gene2coord[record][gene][0])

            gene_median = numpy.median(subset_table.iloc[start_pos:end_pos, 1])
            f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
            record, gene, start_pos, end_pos, gene_median, round(gene_median / median_depth, 2), contig2median_depth[record],
            round(contig2median_depth[record] / median_depth, 2), contig_id2contig_length[record]))

            #print ('all_assembly\t-\t1\t%s\t%s\t-\t-\t-' % (len(samtools_dataframe.iloc[:, 1]),
            #                                                numpy.median(samtools_dataframe.iloc[:, 1])))
