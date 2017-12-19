import json
import mysql.connector
import re
import csv
import itertools

cnx = mysql.connector.connect(option_files=snakemake.params["conf"], option_groups=snakemake.params["id"])
cnx.get_warnings = True
cursor = cnx.cursor()

out = snakemake.output[0]

with open(out, "a") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["Specimen", "software", "gene", "antibiotic_resistance"])
    for j in snakemake.params["soft"]:
        for i in snakemake.params['samples']:
            cmd = "SELECT Specimen, gene from resistance_associated_genes where software=\""+j+"\" and Specimen =\"{0}\";".format(str(i))
            cursor.execute(cmd)
            results = cursor.fetchall()
            for row in results:
                gene = row[1]
                cmd = "SELECT antibiotic from resistance_conferring_genes_annotation where annotation_source=\""+j+"\" and gene = \"{0}\";".format(str(gene))
                cursor.execute(cmd)
                results = cursor.fetchall()
                anti = ",".join([x[0] for x in results])
                res=[i, j, gene, anti]
                w.writerow(res)
                
            

        

cnx.commit()
cnx.close()



