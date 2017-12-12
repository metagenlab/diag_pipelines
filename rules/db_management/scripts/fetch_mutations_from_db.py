import json
import mysql.connector
import re
import csv
import itertools

cnx = mysql.connector.connect(option_files=snakemake.params["conf"], option_groups=snakemake.params["id"])
cnx.get_warnings = True
cursor = cnx.cursor()


for j in snakemake.params["soft"]:
    for i in snakemake.params['samples']:
        cmd = "SELECT * from mutations  where software=\""+j+"\" and Specimen =\"{0}\";".format(str(i))
        cursor.execute(cmd)
        results = cursor.fetchall()
        out = "_".join(snakemake.output[0].split("_")[:-1])+"_"+j+".tsv"
        with open(out, "w") as f:
            print(results)
            w = csv.DictWriter(f, ["Specimen"] + anti, delimiter="\t")
            w.writeheader()
            for key, val in sorted(res.items()):
                row = {"Specimen": key}
                row.update(val)
                w.writerow(row)
cmd="Specimen="
for i in snakemake.params['samples']:
    cmd += '"' + i + '"' + " or Specimen="
cmd = "Select * from mutations  where "+cmd[:-13]+";"

cursor.execute(cmd)
mutations = [x[0] for x in cursor.fetchall()]

with open(snakemake.output[0], "w") as f:
    w = csv.DictWriter(f, ["Specimen"] + anti)
    w.writeheader()
    for key, val in sorted(res.items()):
        row = {"Specimen": key}
        row.update(val)
        w.writerow(row)
    
cnx.commit()
cnx.close()



