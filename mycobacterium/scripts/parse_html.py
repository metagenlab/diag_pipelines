import csv
import os
import mysql.connector
import re
from datetime import datetime
from bs4 import BeautifulSoup as bs
from collections import OrderedDict

cnx = mysql.connector.connect(option_files=snakemake.input[1], option_groups="myco")
cnx.get_warnings = True
cursor = cnx.cursor()
with open(snakemake.output[0], "w") as myfile:
    myfile.write("")

def extract_data(filename, curs, inf, test, anti, log):
    dic={}
    soup=bs(open(filename).read(), "html.parser")
    for mdr in soup.find_all("li"):
        if "Patient: " in mdr.get_text():
            patient="".join(mdr.get_text().replace("Patient:", "").split())
    for lol in soup.find_all("h5"):
        if lol.get_text() =="Genomes":
            table=lol.find_next("table")
            values=table.find_next("tbody").findAll("tr")
            header=table.find_next("thead").findAll("th")
            l=[k.get_text() for k in header]
            for u in values:
                cells=u.findAll("td")
                v=[s.strip() for s in [k.get_text() for k in cells]]
                dic[v[0]]={}
                dic[v[0]]["Patient"]=patient
                for i in range(1,8):
                    if l[i]=="Specimen Collected Date":
                        dic[v[0]]["Specimen_Collected_Date"]=datetime.strptime(v[i], '%b %d, %Y').strftime('%Y-%m-%d')
                    else:
                        dic[v[0]][l[i].replace(" ", "_")]=v[i]
                dic[v[0]]["Filename"]=filename
    if not len(dic):
        print(filename)
        raise ValueError("No genome found in this entry")
    for lol in soup.find_all("h5"):
        if lol.get_text() =="Specimen":
            table=lol.find_next("table")
            header=table.find_next("thead").findAll("th")
            l=[k.get_text() for k in header]
            values=table.find_next("tbody").findAll("tr")
            for u in values:
                cells=u.findAll("td")
                v=[s.strip() for s in [k.get_text() for k in cells]]
                if v[0] in dic.keys():
                    dic[v[0]][l[2]]=v[2]
    for lol in soup.find_all("h5"):
        if lol.get_text() =="Drug susceptibility testing":
            table=lol.find_next("table")
            header=table.find_next("thead").findAll("th")
            l=[k.get_text() for k in header]
            values=table.find_next("tbody").findAll("tr")
            for u in values:
                cells=u.findAll("td")
                v=[s.strip() for s in [k.get_text() for k in cells]]
                if v[0] in dic.keys():
                    dic[v[0]][v[1]]={}
                    for i in range(2, len(v)):
                        if l[i]=="Date":
                            dic[v[0]][v[1]]["Date"]=datetime.strptime(v[i], '%b %d, %Y').strftime('%Y-%m-%d')
                        else:
                            if v[i] != "":
                                dic[v[0]][v[1]][l[i]]=v[i]
    for i in dic.keys():
        cmd=[]
        samp = {k: dic[i][k] for k in inf}
        columns = ', '.join(samp.keys())
        cmd.append("INSERT IGNORE INTO samples ( Specimen, {0} ) VALUES ( \"{1}\", {2} )".format(', '.join(samp.keys()), i, "\""+"\", \"".join(samp.values())+"\""))
        cmd.append("INSERT IGNORE into corres (Specimen, file) VALUES (\"{0}\", \"{1}\")".format(i, re.sub(r'.*/', "", dic[i]["Filename"])))
        if len(dic[i]["Sequence_Read_Archive"].split(",")) > 1:
            for sra in dic[i]["Sequence_Read_Archive"].split(","):
                cmd.append("INSERT IGNORE into sras (Specimen, Sequence_Read_Archive) VALUES (\"{0}\", \"{1}\")".format(i, "".join(sra.split())))
        else:
            cmd.append("INSERT IGNORE into sras (Specimen, Sequence_Read_Archive) VALUES (\"{0}\", \"{1}\")".format(i, dic[i]["Sequence_Read_Archive"]))
        tests_made = [val for val in test if val in list(dic[i].keys())]
        for u in tests_made:
            date = dic[i][u]["Date"]
            for key, value in dic[i][u].items():
                if key != "Date":
                    cmd.append("INSERT IGNORE into dst (Specimen, test, Date, antibio, phenotype) VALUES (\"{0}\", \"{1}\",\"{2}\",\"{3}\", \"{4}\" )".format(i, u, date, key, value))
        for c in cmd:
            curs.execute(c)
            warn=curs.fetchwarnings()
            if warn is not None:
                with open(log, "a") as f:
                    f.write(str(warn[0])+"\n") 
    return(dic)



files = [f for f in os.scandir(snakemake.params[0]) if os.path.isfile(snakemake.params[0]+f.name)]
info = ["Patient", "Bio_Project", "Material", "Specimen_Collected_Date", "Bio_Sample","Lineage","Octal_Spoligotype"]
tests = ["Lowenstein-Jensen", "Bactec", "GeneXpert", "Hain", "DST"]
antibio = ["Specimen","test","Date", 'H', 'R', 'S', 'E', 'Ofx', 'Cm', 'Am', 'Km', 'Z', 'Lfx', 'Mfx', 'Pas', 'Pto', 'Cs', 'Amx/Clv', 'Mb', 'Dld', 'Bdq', 'Ipm/Cln', 'Lzd', 'Cfz', 'Clr', 'Ft', 'AG/CP', 'Action']

for i in files:
    extract_data(snakemake.params[0]+i.name, cursor, info, tests, antibio, snakemake.output[0])

cnx.commit()
cnx.close()
    


