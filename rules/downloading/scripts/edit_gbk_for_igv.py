from Bio import SeqIO 

records = SeqIO.parse(open(snakemake.input[0], 'r'))

out = []
for record in records:
    record.accession = [record.version] 
    out.append(record)
 
SeqIO.write(out, snakemake.output[0], "genbank")
