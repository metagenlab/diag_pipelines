
from Bio import GenBank 

records = GenBank.parse(open(snakemake.input[0], 'r'))

out = []
for record in records:
    record.accession = [record.version] 
    out.append(record)
    
with open(snakemake.output[0], "w") as f:
    for r in out:
        data = str(r).split("\n")
        print(data[-2])
        if 'CONTIG' in data[-2]:
            data[-2] = '//'
            f.write('\n'.join(data[0:-1]) + '\n')
        else:
            f.write('\n'.join(data) + '\n')
