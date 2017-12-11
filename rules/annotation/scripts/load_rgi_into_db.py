import mysql.connector
import csv
import re

cnx = mysql.connector.connect(option_files=snakemake.input[1], option_groups=snakemake.params["db"])
cnx.get_warnings = True
cursor = cnx.cursor()


def load_row_into_db(row, curs, log, sample):
    if row["SNP"] != "n/a":
        if re.search("^\w+ \w+ \w+ mutants conferring resistance to", row["ARO_name"]):
            gene = row["ARO_name"].split()[2]
        elif re.search("^\w+ \w+ \w+ \w+ mutations conferring resistance to", row["ARO_name"]):
            gene = row["ARO_name"].split()[2]+"_"+row["ARO_name"].split()[3]
        elif re.search("^\w+ \w+ \w+ mutations conferring resistance to", row["ARO_name"]):
            gene = row["ARO_name"].split()[2]
        elif re.search("^\w+ \w+ \w+ conferring resistance to", row["ARO_name"]):
            gene = row["ARO_name"].split()[2]
        pos = re.sub("[^0-9]", "", row["SNP"])
        ref = row["SNP"][0]
        mut = row["SNP"][-1]
        try:
            cmd="INSERT IGNORE INTO mutations (Specimen, software, gene, position, ref_aa, mut_aa) VALUES (\"{0}\", \"rgi\", \"{1}\", {2}, \"{3}\", \"{4}\");".format(sample, gene, pos, ref, mut)
            curs.execute(cmd)
            i=curs.fetchwarnings()
            if i is not None:
                with open(log, "a") as f:
                    f.write(str(i)+"\n")

        except UnboundLocalError as err:
            raise ValueError("Failed during the search of the mutated gene in the ARO_name value in the rgi tsv file: {}".format(err))
        for j in row["Best_Hit_ARO_category"].strip().split(";"):
            if "resistance" in j:
                if re.search("determinant of \w+(?:\S\w+) resistance", j):
                    antibio = j.split(" ")[-2]
                elif re.search("determinant of resistance to \w+(?:\S\w+) antibiotics", j):
                    antibio = j.split(" ")[-2]+"_antibiotics"
                try:
                    cmd="INSERT IGNORE INTO phenotype_prediction (Specimen, software, antibiotic) VALUES (\"{0}\", \"rgi\", \"{1}\");".format(sample, antibio)
                    curs.execute(cmd)
                    i=curs.fetchwarnings()
                    if i is not None:
                        with open(log, "a") as f:
                            f.write(str(i)+"\n")
                except UnboundLocalError as err:
                    raise ValueError("Failed during the search of the antibiotic resistance in the Best_Hit_ARO_category value, for a resistance conferring gene mutation, in the rgi tsv file: {}".format(err))

                
        
    elif "intrinsic" not in row["Best_Hit_ARO"]:
        gene = row["Best_Hit_ARO"].replace(" ", "_")
        cmd="INSERT IGNORE INTO resistance_conferring_genes (Specimen, software, gene, description) VALUES (\"{0}\", \"rgi\", \"{1}\", \"{2}\");".format(sample, gene, row["Best_Hit_ARO_category"].strip())
        curs.execute(cmd)
        i=curs.fetchwarnings()
        if i is not None:
            with open(log, "a") as f:
                f.write(str(i)+"\n")
        for j in row["Best_Hit_ARO_category"].strip().split(";"):
            if "resistance" in j and "efflux pump" not in j and "gene cluster" not in j:
                if re.search("determinant of \w+(?:\S\w+) resistance", j):
                    antibio = j.split(" ")[-2]
                elif re.search("determinant of resistance to \w+(?:\S\w+) antibiotics", j):
                    antibio = j.split(" ")[-2]+"_antibiotics"
                try:
                    cmd="INSERT IGNORE INTO phenotype_prediction (Specimen, software, antibiotic) VALUES (\"{0}\", \"rgi\", \"{1}\");".format(sample, antibio)
                    curs.execute(cmd)
                    i=curs.fetchwarnings()
                    if i is not None:
                        with open(log, "a") as f:
                            f.write(str(i)+"\n")
                except UnboundLocalError as err:
                    raise ValueError("Failed during the search of the antibiotic resistance in the Best_Hit_ARO_category value, for a resistance conferring gene acquisition, in the rgi tsv file: {}".format(err))

                
            

with open(snakemake.output[0], "w") as myfile:
    myfile.write("")

with open(snakemake.input[0], "r") as tsvfile:
    csvreader = csv.DictReader(tsvfile, delimiter="\t")
    for row in csvreader:
        load_row_into_db(row, cursor, snakemake.output[0], snakemake.params["id"][snakemake.wildcards["sample"]])
        
cnx.commit()
cnx.close()
