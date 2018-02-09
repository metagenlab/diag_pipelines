import pandas
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna, generic_protein

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]

snps = pandas.read_csv(snakemake.input["genotype"], sep="\t", header=0)
snps.rename(axis="columns", mapper= lambda x: re.sub("\[[0-9]+\]", "", x.replace(":GT", "")), inplace=True)
snps.replace(to_replace=".", value=0, inplace=True)

print(snps.columns)

ref = snakemake.wildcards["ref"]
all_samples = snakemake.params["samples"]
all_samples.sort(key=natural_keys)

ref_sequence = list(snps.loc[:, "REF"])
print(len(ref_sequence))

res={}
for sample in all_samples:
    res[sample] = list(ref_sequence)

position = 0
for index, row in snps.iterrows():
    alts = row["ALT"].split(",")
    for sample in all_samples:
        if int(row[sample]):
            res[sample][position] = alts[int(row[sample])-1]
    position += 1



with open(snakemake.output[0], "a") as align_file:
    for sample in res.keys():
        record = SeqRecord(Seq("".join(res[sample]), generic_dna), id=sample, description="")
        SeqIO.write(record, align_file, "fasta")
    record = SeqRecord(Seq("".join(ref_sequence), generic_dna), id=ref, description="")
    SeqIO.write(record, align_file, "fasta")
    
