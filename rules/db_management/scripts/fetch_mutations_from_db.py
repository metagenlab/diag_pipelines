import json
import mysql.connector
import re
import csv
import itertools
import pandas

cnx = mysql.connector.connect(option_files=snakemake.params["conf"], option_groups=snakemake.params["id"])
cnx.get_warnings = True
cursor = cnx.cursor()

out_csv = snakemake.output[0]
out_xlsx = snakemake.output[1]

spec = snakemake.params["sample_corres"][snakemake.wildcards["sample"]]



all_res = []

with open(out_csv, "a") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["Software", "Gene", "Resistance Type", "Variant", "Antibiotic resistance prediction"])
    for j in snakemake.params["soft"]:
        cmd = "SELECT * from resistance_associated_mutations where software=\""+j+"\" and Specimen =\""+spec+"\";"
        cursor.execute(cmd)
        results = cursor.fetchall()
        for result in results:
            gene = result[2]
            cmd = "SELECT distinct antibiotic from resistance_conferring_mutations_annotation where annotation_source=\""+j+"\" and gene=\""+gene+"\";"
            cursor.execute(cmd)
            annot = cursor.fetchall()
            if len(annot)>1:
                raise ValueError("One gene is annotated for more than one antibiotic in the mutation db")
            else:
                anti = annot[0][0]
                mut = result[4]+str(result[3])+result[5]
                res = [j, gene, "variant", mut, anti]
                w.writerow(res)
                all_res.append(res)
        cmd = "SELECT Specimen, gene from resistance_associated_genes where software=\""+j+"\" and Specimen =\""+spec+"\";"
        cursor.execute(cmd)
        results = cursor.fetchall()
        for row in results:
            gene = row[1]
            cmd = "SELECT antibiotic from resistance_conferring_genes_annotation where annotation_source=\""+j+"\" and gene = \""+str(gene)+"\";"
            cursor.execute(cmd)
            results = cursor.fetchall()
            anti = ",".join([x[0] for x in results])
            res=[j, gene, "gene presence", "",  anti]
            w.writerow(res)
            all_res.append(res)
                
cnx.commit()
cnx.close()

df = pandas.DataFrame(all_res, columns= ["Software", "Gene", "Resistance Type", "Variant", "Antibiotic resistance prediction"])


writer = pandas.ExcelWriter(snakemake.output[1])
df.to_excel(writer, snakemake.params["sample_corres"][snakemake.wildcards["sample"]], index=False)
writer.save()
