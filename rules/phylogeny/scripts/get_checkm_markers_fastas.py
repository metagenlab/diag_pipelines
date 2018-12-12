
import pandas
import re
from Bio import SeqIO

checkm_table = snakemake.input["checkm_table"]
checkm_fastas = snakemake.input["checkm_fastas"]


def parse_checkm_marker_table(checkm_table):
    marker_table = pandas.read_csv(checkm_table, delimiter='\t', header=0, index_col=0)

    marker2sample2records = {}

    for n, row in marker_table.iterrows():
        sample_id = row.name
        marker = row[0]
        gene_id = row[1]
        if marker not in marker2sample2records:
            marker2sample2records[marker] = {}
        if sample_id not in marker2sample2records[marker]:
            marker2sample2records[marker][sample_id] = [gene_id]
        else:
            marker2sample2records[marker][sample_id].append(gene_id)
    # remove hits with multiple copies
    for marker in marker2sample2records:
        for sample in marker2sample2records[marker]:
            if len(marker2sample2records[marker][sample]) > 1:
                marker2sample2records[marker].pop(sample)
    return marker2sample2records




def parse_fastas(fasta_files):

    sample2genes = {}

    for fasta_file in fasta_files:
        sample = fasta_file.split("/")[4]
        sample_dico = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
        sample2genes[sample] = sample_dico
    return sample2genes

marker_genes = parse_checkm_marker_table(checkm_table)
fasta_dico = parse_fastas(checkm_fastas)

marker2records = {}
for marker in marker_genes:
    marker_records = []
    for sample in marker_genes[marker]:
        print("%s\t%s\t%s" % (marker, sample, marker_genes[marker][sample]))
        record = fasta_dico[str(sample)][marker_genes[marker][sample][0]]
        record.name = str(sample)
        record.id = str(sample)
        record.description = ""
        record.seq = record.seq[0:-1]
        marker_records.append(record)
    with open("phylogeny/checkm/marker_fastas/%s.faa" % marker, "w") as f:
        print(marker_records)
        SeqIO.write(marker_records, f, "fasta")
