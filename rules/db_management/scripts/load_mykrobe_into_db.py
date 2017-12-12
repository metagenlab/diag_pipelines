import csv
import mysql.connector

cnx = mysql.connector.connect(option_files=snakemake.params["conf"], option_groups=snakemake.params["db"])
cnx.get_warnings = True
cursor = cnx.cursor()


def load_row_into_db(row, curs, log, sample):
    cmds = []
    if row["susceptibility"] == "R":
        gene = row["genes (prot_mut-ref_mut:percent_covg:depth)"].split(":")[0]
        anti = row["drug"].lower()
        cmds.append("INSERT INTO antibiotic_resistance_conferring_genes_annotation (gene, annotation_source, antibiotic) VALUES (\"{0}\", \"mykrobe\", \"{1}\");".format(gene, anti))
        cmds.append("INSERT INTO phenotype_prediction (specimen, software, antibiotic) VALUES (\"{0}\", \"mykrobe\", \"{1}\");".format(sample, anti))
        cmds.append("INSERT INTO resistance_associated_genes(specimen, software, gene) VALUES  (\"{0}\", \"mykrobe\", \"{1}\");".format(sample, gene))

        for cmd in cmds:
            try:
                cursor.execute(cmd)
            except mysql.connector.errors.Error as err :
                with open(log, "a") as f:
                        f.write("Something went wrong: {}\n".format(err))




with open(snakemake.output[0], "w") as logfile:
    logfile.write("")

with open(snakemake.input[0], "r") as tsvfile:
    csvreader = csv.DictReader(tsvfile, delimiter="\t")
    for row in csvreader:
        load_row_into_db(row, cursor, snakemake.output[0], snakemake.params["id"][snakemake.wildcards["sample"]])



cnx.commit()
cnx.close()
