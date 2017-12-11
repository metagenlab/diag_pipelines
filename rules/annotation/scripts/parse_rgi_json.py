import json
import mysql.connector
import pronto

cnx = mysql.connector.connect(option_files=snakemake.input[1], option_groups=snakemake.params["db"])
cnx.get_warnings = True
cursor = cnx.cursor()

print(ont.obo)

def parse_rgi_json(filename, curs, log, sample):
    with open(filename, "r") as f:
        rgi = json.load(f)
        for j in rgi.keys():
            if j != "_metadata":
                for k in rgi[j].keys():
                    aro_accession = rgi[j][k]["ARO_acession"]
                    
                    try: 
                        snp=rgi[j][k]["SNP"]
                        try:
                            antibio=rgi[j][k]["ARO_name"].split()[-1].strip()
                            gene = rgi[j][k]["ARO_name"].split()[2]
                            pos = snp["position"]
                            ref = snp["original"]
                            mut = snp["change"]
                            cmd="INSERT IGNORE INTO mutations (Specimen, software, gene, position, ref_aa, mut_aa) VALUES (\"{0}\", \"rgi\", \"{1}\", {2}, \"{3}\", \"{4}\");".format(str(sample), gene, pos, ref, mut)
                            curs.execute(cmd)
                            i=curs.fetchwarnings()
                            if i is not None:
                                with open(log, "a") as f:
                                    f.write(i[2]+"\n")
                            cmd="INSERT IGNORE INTO phenotype_prediction (Specimen, software, antibiotic) VALUES (\"{0}\", \"rgi\", \"{1}\");".format(str(sample), str(antibio))
                            curs.execute(cmd)
                            i=curs.fetchwarnings()
                            if i is not None:
                                with open(log, "a") as f:
                                    f.write(str(i)+"\n")
                    
                        except IndexError:
                            print(rgi[j][k]["SNP"])
                            pass
                    except KeyError:
                        if rgi[j][k]["model_type"] == "protein homolog model":
                            gene=rgi[j][k]["ARO_name"].replace(" ", "_")
                            cmd="INSERT IGNORE INTO resistance_confering_genes (Specimen, software, gene) VALUES (\"{0}\", \"rgi\", \"{1}\");".format(str(sample), str(gene))
                            curs.execute(cmd)
                            i=curs.fetchwarnings()
                            if i is not None:
                                with open(log, "a") as f:
                                    f.write(str(i)+"\n")
                            for i in rgi[j][k]["ARO_category"].keys():
                                print(rgi[j][k]["ARO_category"][i]["category_aro_name"])


                            
    return("ok")


with open(snakemake.output[0], "w") as myfile:
    myfile.write("")

parse_rgi_json(snakemake.input[0], cursor, snakemake.output[0], snakemake.wildcards.sample)
cnx.commit()
cnx.close()
