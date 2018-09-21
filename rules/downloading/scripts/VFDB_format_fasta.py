#!/usr/bin/python
from Bio import SeqIO


fasta_file = snakemake.input[0]
out_fasta = snakemake.output[0]
out_table = snakemake.output[1]

f = open(out_fasta, 'w')
t = open(out_table, 'w')

reformated_records = []


for entry in SeqIO.parse(open(fasta_file, "rU"), "fasta"):
    try:
        seq_accession = entry.name.split('|')[1][0:-1]
        seq_id = entry.name.split('|')[0][0:-3]
    except IndexError:
        seq_accession = '-'
        seq_id = entry.name

    description = entry.description

    try:
        gene = description.split(') (')[1].split(' ')[0][0:-1]
        product = description.split(')')[2].split('[')[0][1:]
    except IndexError:
        gene = description.split('(')[1].split(' ')[0][0:-1]
        product = description.split(')')[1].split('[')[0][1:]
    vf_description = description.split('[')[1][0:-2]
    vf_name = vf_description.split(' (')[0]

    # attention VFG000881 pas d'accession
    # edit du fichier fasta => blast identique: WP_000930062.1
    try:
        gbk_accession = entry.name.split('gb|')[1].split(')')[0]
    except:
        gbk_accession = '-'
    # handle the 2 different fasta headers
    # 1. >VFG007172(gb|NP_798073) (vscF) type III secretion system needle protein VscF [T3SS1 (VF0408)] [Vibrio parahaemolyticus RIMD 2210633]
    # 2.
    if description.count('[') == 3:
        vf_description = description.split(' [')[1][0:-2]
        vf_id = vf_description.split(' (')[1][0:]
        bacteria = description.split(' [')[2][0:-1]
        product = description.split(') ')[2].split('] [')[0].split(' [')[0]
    else:
        vf_description = description.split('[')[1][0:-2]
        vf_id = vf_description.split(' (')[1][0:-1]
        bacteria = description.split('[')[2][0:-1]

    entry.id = seq_id
    entry.description = vf_description
    reformated_records.append(entry)
    t.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (seq_id, vf_id, gene, bacteria, vf_description.split('(')[0], product, gbk_accession))

SeqIO.write(reformated_records, f, "fasta")
t.close()
f.close()
