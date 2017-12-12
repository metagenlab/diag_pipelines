import csv
import mysql.connector

cnx = mysql.connector.connect(option_files=snakemake.params["conf"], option_groups=snakemake.params["db"])
cnx.get_warnings = True
cursor = cnx.cursor()


def load_row_into_db(row, curs, log, sample):
    gene = row[4]
    cmd="INSERT IGNORE INTO resistance_associated_genes(specimen, software, gene) VALUES (\"{0}\", \"abricate\", \"{1}\");".format(sample, str(gene))
    curs.execute(cmd)
    i=cursor.fetchwarnings()
    if i is not None:
        with open(log, "a") as f:
            f.write(str(i)+"\n")

with open(snakemake.output[0], "w") as logfile:
    logfile.write("")


with open(snakemake.input[0]) as abricate_file:
    csvreader = csv.reader(abricate_file, delimiter='\t')
    next(csvreader, None)
    for row in csvreader:
        load_row_into_db(row, cursor, snakemake.output[0], snakemake.params["id"][snakemake.wildcards["sample"]])
    


cnx.commit()
cnx.close()
