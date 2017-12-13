import csv
import mysql.connector
import re

cnx = mysql.connector.connect(option_files=snakemake.params["conf"], option_groups=snakemake.params["db"])
cnx.get_warnings = True
cursor = cnx.cursor()


def load_row_into_db(row, curs, log, sample):
    cmds = []
    if row["susceptibility"] == "R" and row["variants (gene:alt_depth:wt_depth:conf)"] =="":
        gene = row["genes (prot_mut-ref_mut:percent_covg:depth)"].split(":")[0]
        anti = row["drug"].lower()
        cmds.append("INSERT INTO resistance_conferring_genes_annotation (gene, annotation_source, antibiotic) VALUES (\"{0}\", \"mykrobe\", \"{1}\");".format(gene, anti))
        cmds.append("INSERT INTO phenotype_prediction_from_gene_presence (specimen, software, antibiotic) VALUES (\"{0}\", \"mykrobe\", \"{1}\");".format(sample, anti))
        cmds.append("INSERT INTO resistance_associated_genes(specimen, software, gene) VALUES  (\"{0}\", \"mykrobe\", \"{1}\");".format(sample, gene))
        
    if row["susceptibility"] == "R" and row["variants (gene:alt_depth:wt_depth:conf)"] !="" and "del" not in row["variants (gene:alt_depth:wt_depth:conf)"] and "ins" not in row["variants (gene:alt_depth:wt_depth:conf)"]:
        gene = row["variants (gene:alt_depth:wt_depth:conf)"].split('_')[0]
        anti = row["drug"].lower()
        mutation = row["variants (gene:alt_depth:wt_depth:conf)"].split('_')[1].split("-")[0]
        ref = mutation[0]
        mut = mutation[-1]
        pos = re.sub("[^0-9]", "", mutation)
        cmds.append("INSERT INTO resistance_associated_mutations (specimen, software, gene, position, ref_aa, mut_aa) VALUES (\"{0}\", \"mykrobe\", \"{1}\", {2}, \"{3}\", \"{4}\");".format(sample, gene, pos, ref, mut))
        cmds.append("INSERT INTO phenotype_prediction_from_mutation (specimen, software, antibiotic) VALUES (\"{0}\", \"mykrobe\", \"{1}\");".format(sample, anti))
        cmds.append("INSERT INTO resistance_conferring_mutations_annotation (gene, position, ref_aa, mut_aa, annotation_source, antibiotic) VALUES (\"{0}\", {1}, \"{2}\", \"{3}\", \"mykrobe\", \"{4}\");".format(gene, pos, ref, mut, anti))
    if "del" in row["variants (gene:alt_depth:wt_depth:conf)"] or "ins" in row["variants (gene:alt_depth:wt_depth:conf)"]:
        raise Exception("Insertion or deletion confering resistance have been detected by Mykrobe in {0}, parsing currently not implemented".format(sample))
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
