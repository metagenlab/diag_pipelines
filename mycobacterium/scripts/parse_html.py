from bs4 import BeautifulSoup as bs
import csv
from collections import OrderedDict
from dateutil import parser
import mysql.connector

def extract_data(filename):
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
                        dic[v[0]]["Specimen_Collected_Date"]=parser.parse(v[i]).strftime('%Y-%m-%d')
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
                            dic[v[0]][v[1]]["Date"]=parser.parse(v[i]).strftime('%Y-%m-%d')
                        else:
                            dic[v[0]][v[1]][l[i]]=v[i]

    return(dic)

final_dic={}
files = [f for f in os.scandir("../") if os.path.isfile("../"+f.name)]
print(files)
for i in files:
    final_dic={** final_dic, ** extract_data("../"+i.name)}
    
     

tests = ["Lowenstein-Jensen", "Bactec", "GeneXpert", "Hain", "DST"]
info = ["Specimen","Patient","Bio_Project","Material","Specimen_Collected_Date", "Bio_Sample","Lineage","Octal_Spoligotype"]
antibio = ["Specimen","test","Date", 'H', 'R', 'S', 'E', 'Ofx', 'Cm', 'Am', 'Km', 'Z', 'Lfx', 'Mfx', 'Pas', 'Pto', 'Cs', 'Amx/Clv', 'Mb', 'Dld', 'Bdq', 'Ipm/Cln', 'Lzd', 'Cfz', 'Clr', 'Ft', 'AG/CP', 'Action']

deli="\t"
with open("samples.csv", "w") as f:
    w = csv.DictWriter(f, None, delimiter=deli)
    w.fieldnames = info
    w.writeheader()


with open("test.csv", "w") as f:
    w = csv.DictWriter(f, None, delimiter=deli)
    w.fieldnames = antibio
    w.writeheader()

with open("sras.csv", "w") as f:
    w = csv.DictWriter(f, None, delimiter=deli)
    w.fieldnames=["Specimen", "Sequence_Read_Archive"]
    w.writeheader()

print(final_dic.keys())
with open('samples.csv', 'a') as f:
    for i in final_dic.keys():
        tmp = {}
        tmp["Specimen"]=i
        tmp = {** tmp, ** {k: final_dic[i][k] for k in info[1:]}}
        w = csv.DictWriter(f, tmp, delimiter=deli)
        w.fieldnames=info
        w.writerow(tmp)
        with open("sras.csv", "a") as u:
            for sras in final_dic[i]["Sequence_Read_Archive"].split(","):
                u.write(i+deli+sras.replace(" ", "")+"\n")
        tests_made = [val for val in tests if val in list(final_dic[i].keys())]
        with open("test.csv", "a") as g:
            for test in tests_made:
                tmpres = {}
                tmpres["Specimen"] = i
                tmpres["test"] = test
                tmpres = {** tmpres, ** final_dic[i][test]}
                print(tmpres)
                w = csv.DictWriter(g, tmpres, delimiter=deli)
                w.fieldnames=antibio
                w.writerow(tmpres)
        
        
        
with open("file_correspondance.csv", "w") as f:
    f.write("Speciment\tID\n")
    for i in final_dic.keys():
        f.write(i+deli+final_dic[i]["Filename"].replace(".", "").replace("/", "")+"\n")

