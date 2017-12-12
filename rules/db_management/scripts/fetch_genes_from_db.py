import json
import mysql.connector
import re
import csv
import itertools

cnx = mysql.connector.connect(option_files=snakemake.params["conf"], option_groups=snakemake.params["id"])
cnx.get_warnings = True
cursor = cnx.cursor()

for j in snakemake.params["soft"]:
    out = "_".join(snakemake.output[0].split("_")[:-1])+"_"+j+".tsv"
    with open(out, "w") as f:
        w = csv.writer(f)
        w.writerow(["Specimen", "gene"])
        for i in snakemake.params['samples']:
            cmd = "SELECT Specimen, gene from resistance_conferring_genes where software=\""+j+"\" and Specimen =\"{0}\";".format(str(i))
            cursor.execute(cmd)
            results = cursor.fetchall()
            for row in results:
                w.writerow(row)
                
            

        

cnx.commit()
cnx.close()



