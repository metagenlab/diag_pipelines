import json
import mysql.connector
import re
import csv
import itertools

cnx = mysql.connector.connect(option_files=snakemake.input[0], option_groups="myco")
cnx.get_warnings = True
cursor = cnx.cursor()

res = {}

cmd="Specimen="
for i in snakemake.params['samples']:
    cmd += '"' + i + '"' + " or Specimen="

cmd = "Select distinct antibiotic from phenotype_prediction where "+cmd[:-13]+";"

cursor.execute(cmd)
anti = [x[0] for x in cursor.fetchall()]


for i in snakemake.params['samples']:
    cmd = "SELECT antibiotic from phenotype_prediction where Specimen =\"{0}\";".format(str(i))
    cursor.execute(cmd)
    res[i]={}
    for k in anti:
        res[i][k]="S"
    for k in [x[0] for x in cursor.fetchall()]:
        res[i][k]="R"


with open(snakemake.output[0], "w") as f:
    w = csv.DictWriter(f, ["Specimen"] + anti)
    w.writeheader()
    for key, val in sorted(res.items()):
        row = {"Specimen": key}
        row.update(val)
        w.writerow(row)
    
cnx.commit()
cnx.close()



