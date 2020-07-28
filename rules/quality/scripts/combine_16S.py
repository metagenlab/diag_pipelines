import re
from Bio import SeqIO

input_files = snakemake.input

rg_node_name = 'NODE_[0-9]+'
rg_node_length= 'length_([0-9]+)'
rg_cov= 'cov_([0-9]+)'

edited_records = []
for one_file in input_files:
    sample_name = one_file.split("/")[0]
    records = SeqIO.parse(one_file, "fasta")
    for record in records:
        # >16S_rRNA::NODE_1353_length_1702_cov_240.973968:15-1553(-)
        node_name = re.search(rg_node_name, record.description).group(0)
        node_length = re.search(rg_node_length, record.description).group(1)
        cov = re.search(rg_cov, record.description).group(1)
        record.name = ''
        record.description = ''
        record.id = f'{sample_name}|{node_name}|length_{node_length}|cov_{cov}'
        edited_records.append(record)

SeqIO.write(edited_records, snakemake.output[0], "fasta")
