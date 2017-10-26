from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import mysql.connector

mysql_option_files = [
    '/home/sacha/.my.cnf',
]

cnx = mysql.connector.connect(option_files=mysql_option_files, option_groups="myco")
cursor = cnx.cursor()

cursor.execute("SELECT gene FROM resistance;")

genes = set([])
for i in cursor:
    genes.add(i[0])
cursor.close()
cnx.close()

print(len(genes))

table = 11

def extract_genes(filename, tot):
    cds = []
    rna = []
    found = set([])
    print(filename)
    for record in SeqIO.parse(filename, "genbank"):
        for f in record.features:
            if f.type=="rRNA":
                if f.qualifiers["locus_tag"][0] in tot:
                    found.add(f.qualifiers["locus_tag"][0])
                    rna.append(f)
                else:
                    try:
                        if f.qualifiers["gene"][0] in tot:
                            found.add(f.qualifiers["gene"][0])
                            rna.append(f)
                    except KeyError:
                        pass

            if f.type=="CDS":
                if f.qualifiers["locus_tag"][0] in tot:
                    found.add(f.qualifiers["locus_tag"][0])
                    cds.append(f)
                    
                else:
                    try:
                        if f.qualifiers["gene"][0] in tot:
                            found.add(f.qualifiers["gene"][0])
                            cds.append(f)
                    except KeyError:
                        pass
    return(found, cds, rna)

def get_prot_seq(record, table=11):
    try:
        return(SeqRecord(id=record.qualifiers["gene"][0], description="H37RV", seq=Seq(record.qualifiers["translation"][0])))
    except KeyError:
        return(SeqRecord(id=record.qualifiers["locus_tag"][0], description="H37RV", seq=Seq(record.qualifiers["translation"][0])))

            

def get_nucl_seq(record, gb_record):
    try:
        return(SeqRecord(id=record.qualifiers["gene"][0], seq=record.extract(gb_record.seq), description="H37RV"))
    except KeyError:
        return(SeqRecord(id=record.qualifiers["locus_tag"][0], seq=record.extract(gb_record.seq), description="H37RV"))

        
k = extract_genes(snakemake.input[0], genes)#, snakemake.output[0], snakemake.output[1])
genome = SeqIO.read(snakemake.input[0], "genbank")

SeqIO.write((get_prot_seq(r) for r in k[1]), snakemake.output[1], "fasta")
SeqIO.write((get_nucl_seq(r, genome) for r in k[1]+k[2]), snakemake.output[0], "fasta")

print(genes-k[0])
