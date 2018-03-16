import mysql.connector
import csv
import re
import pronto

cnx = mysql.connector.connect(option_files=snakemake.params["conf"], option_groups=snakemake.params["db"])
cnx.get_warnings = True
cursor = cnx.cursor()

aro = pronto.Ontology(snakemake.params["ontology_aro"])
ro = pronto.Ontology(snakemake.params["ontology_ro"])
mo = pronto.Ontology(snakemake.params["ontology_mo"])
currated_genes = snakemake.params["currated_genes"]

aro.merge(ro)
aro.merge(mo)


def load_row_into_db(row, curs, log, sample, aro_ont, gene_list):
    cmds = []
    if row["SNP"] != "n/a":
        best_hit = row["Best_Hit_ARO"].split()
        genes = list(set(best_hit).intersection(set(gene_list)))
        if len(genes) == 1:
            gene = genes[0]
            terms = row["ARO"].split(",")
            SNPS = row["SNP"].split(",")
            for i in range(len(terms)):
                term = aro_ont[terms[i].strip()]
                snp = SNPS[i].strip()
                pos = re.sub("[^0-9]", "", snp)
                ref = snp[0]
                mut = snp[-1]
                for child in term.relations:
                    if child.obo_name=="confers_resistance_to_drug" or  child.obo_name=="confers_resistance_to":
                        for antibiotic in term.relations[child]:
                            anti = antibiotic.name.replace("antibiotic", "").strip()
                            cmds.append("INSERT IGNORE INTO resistance_associated_mutations (specimen, software, gene, position, ref_aa, mut_aa) VALUES (\"{0}\", \"rgi\", \"{1}\", {2}, \"{3}\", \"{4}\");".format(sample, gene, pos, ref, mut))
                            cmds.append("INSERT IGNORE INTO phenotype_prediction_from_mutation (specimen, software, antibiotic) VALUES (\"{0}\", \"rgi\", \"{1}\");".format(sample, anti))
                            cmds.append("INSERT IGNORE INTO resistance_conferring_mutations_annotation (gene, position, ref_aa, mut_aa, annotation_source, antibiotic) VALUES (\"{0}\", {1}, \"{2}\", \"{3}\", \"rgi\", \"{4}\");".format(gene, pos, ref, mut, anti))
            for cmd in cmds:
                try:
                    cursor.execute(cmd)
                except mysql.connector.errors.Error as err :
                    with open(log, "a") as f:
                        f.write("Something went wrong: {}\n".format(err))
        else:
            raise Exception("Problem parsing the predicted ARO Model from RGI from {0}".format(sample))
    elif "intrinsic" not in row["Best_Hit_ARO"]:
        term = aro_ont[row["ARO"]]
        gene = row["Best_Hit_ARO"].replace(" ", "_").strip()
        print(term)
        for child in term.relations:
            if child.obo_name=="confers_resistance_to_drug" or  child.obo_name=="confers_resistance_to":
                for antibiotic in term.relations[child]:
                    anti = antibiotic.name.replace("antibiotic", "").strip()
                    cmds.append("INSERT IGNORE INTO phenotype_prediction_from_gene_presence (specimen, software, antibiotic) VALUES (\"{0}\", \"rgi\", \"{1}\");".format(sample, anti))
                    cmds.append("INSERT IGNORE INTO resistance_conferring_genes_annotation (gene, annotation_source, antibiotic) VALUES (\"{0}\", \"rgi\", \"{1}\");".format(gene, anti))
        
        for parent in term.rparents():
            print(parent)
            for child in parent.relations:
                if child.obo_name=="confers_resistance_to_drug" or  child.obo_name=="confers_resistance_to":
                    for antibiotic in parent.relations[child]:
                        anti = antibiotic.name.replace("antibiotic", "").strip()
                        print(anti)
                        cmds.append("INSERT IGNORE INTO phenotype_prediction_from_gene_presence (specimen, software, antibiotic) VALUES (\"{0}\", \"rgi\", \"{1}\");".format(sample, anti))
                        cmds.append("INSERT IGNORE INTO resistance_conferring_genes_annotation (gene, annotation_source, antibiotic) VALUES (\"{0}\", \"rgi\", \"{1}\");".format(gene, anti))

        cmds.append("INSERT IGNORE INTO resistance_associated_genes(specimen, software, gene) VALUES  (\"{0}\", \"rgi\", \"{1}\");".format(sample, gene))
        for cmd in cmds:
            try:
                cursor.execute(cmd)
            except mysql.connector.errors.Error as err :
                with open(log, "a") as f:
                    f.write("Something went wrong: {}\n".format(err))
                    
with open(snakemake.output[0], "w") as myfile:
    myfile.write("")

with open(snakemake.input[0], "r") as tsvfile:
    csvreader = csv.DictReader(tsvfile, delimiter="\t")
    for row in csvreader:
        load_row_into_db(row, cursor, snakemake.output[0], snakemake.params["id"][snakemake.wildcards["sample"]], aro, currated_genes)
        
cnx.commit()
cnx.close()
