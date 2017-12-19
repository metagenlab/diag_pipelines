import json
import mysql.connector
import re
import csv
import itertools

cnx = mysql.connector.connect(option_files=snakemake.params["conf"], option_groups=snakemake.params["id"])
cnx.get_warnings = True
cursor = cnx.cursor()

res = {}

for j in snakemake.params["soft"]:
    cmd="Specimen="
    for i in snakemake.params['samples']:
        cmd += '"' + i + '"' + " or Specimen="
    cmd = "Select distinct antibiotic from phenotype_prediction_from_gene_presence where software=\""+j+"\" and ("+cmd[:-13]+");"
    cursor.execute(cmd)
    anti = [x[0] for x in cursor.fetchall()]
    for i in snakemake.params['samples']:
        cmd = "SELECT antibiotic from phenotype_prediction_from_gene_presence  where software=\""+j+"\" and Specimen =\"{0}\";".format(str(i))
        cursor.execute(cmd)
        res[i]={}
        for k in anti:
            res[i][k]="S"
        for k in [x[0] for x in cursor.fetchall()]:
            res[i][k]="R"
    out = "_".join(snakemake.output[0].split("_")[:-1])+"_"+j+".tsv"
    with open(out, "w") as f:
        w = csv.DictWriter(f, ["Specimen"] + anti, delimiter="\t")
        w.writeheader()
        for key, val in sorted(res.items()):
            row = {"Specimen": key}
            row.update(val)
            w.writerow(row)

            
cnx.commit()
cnx.close()



