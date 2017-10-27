import vcf
import csv
from Bio import SeqIO
from Bio.Seq import Seq
import mysql.connector

cnx = mysql.connector.connect(option_files=snakemake.input[2], option_groups="myco")
cnx.get_warnings = True
cursor = cnx.cursor()


vcf_reader = vcf.Reader(open(snakemake.input[0], 'r'))
genes ={}

for index, record in enumerate(SeqIO.parse(snakemake.input[1], "fasta")):
    genes[record.id]=record

with open(snakemake.output[0], "w") as myfile:
    myfile.write("")
    
for record in vcf_reader:
    mutation = str(record.ALT).replace("[","").replace("]","")
    print(record)
    if mutation in ["A","T","G","C"]:
        codon = genes[record.CHROM][((record.POS-1)//3*3):(((record.POS-1)//3+1)*3)].seq
        mutated_codon=list(codon)
        mutated_codon[record.POS%3-1]=mutation
        mutated_codon=''.join(mutated_codon)
        if str(codon.translate()) != str(Seq(mutated_codon).translate()):
            cmd="INSERT IGNORE INTO mutations (Specimen, software, gene, position, ref_aa, mut_aa) VALUES (\"{0}\", \"local\", \"{1}\", {2}, \"{3}\", \"{4}\");".format(snakemake.wildcards.sample, record.CHROM, record.POS//3, str(codon.translate()), str(Seq(mutated_codon).translate()))
            cursor.execute(cmd)
            print(cmd)
            i=cursor.fetchwarnings()
            if i is not None:
                with open(snakemake.output[0], "a") as f:
                    f.write(str(i)+"\n")
    else:
        print("INDEL, do something else!")
    
cnx.commit()
cnx.close()
