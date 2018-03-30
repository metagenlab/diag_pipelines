import pandas
import pronto
import csv
import re

aro = pronto.Ontology(snakemake.params["ontology_aro"])
ro = pronto.Ontology(snakemake.params["ontology_ro"])
mo = pronto.Ontology(snakemake.params["ontology_mo"])


aro.merge(ro)
aro.merge(mo)


def extract_rgi(row, aro_ont, gene_list):
    res = []
    if row["SNP"] != "n/a":
        pass
#        best_hit = row["Best_Hit_ARO"].split()
#        genes = list(set(best_hit).intersection(set(gene_list)))
#        if len(genes) == 1:
#            gene = genes[0]
#            terms = row["ARO"].split(",")
#            SNPS = row["SNP"].split(",")
#            for i in range(len(terms)):
#                term = aro_ont[terms[i].strip()]
#                snp = SNPS[i].strip()
#                pos = re.sub("[^0-9]", "", snp)
#                ref = snp[0]
#                mut = snp[-1]
#                for child in term.relations:
#                    if child.obo_name=="confers_resistance_to_drug" or  child.obo_name=="confers_resistance_to":
#                        for antibiotic in term.relations[child]:
#                            anti = antibiotic.name.replace("antibiotic", "").strip()
#                            res.append(["rgi", gene, "variant", ref+str(pos)+mut, anti])
#        else:
#            raise Exception("Problem parsing the predicted ARO Model from RGI")
    elif "intrinsic" not in row["Best_Hit_ARO"]:
        for terms in row["ARO"].split(","):
            term = aro_ont[terms.strip()]
            gene = row["Best_Hit_ARO"].replace(" ", "_").strip()
            for child in term.relations:
                if child.obo_name=="confers_resistance_to_drug" or  child.obo_name=="confers_resistance_to":
                    for antibiotic in term.relations[child]:
                        anti = antibiotic.name.replace("antibiotic", "").strip()
                        res.append(["rgi", gene, "gene presence", "", anti])
            for parent in term.rparents():
                for child in parent.relations:
                    if child.obo_name=="confers_resistance_to_drug" or  child.obo_name=="confers_resistance_to":
                        for antibiotic in parent.relations[child]:
                            anti = antibiotic.name.replace("antibiotic", "").strip()
                            res.append(["rgi", gene, "gene presence", "", anti])
    return(res)

all_res = []


with open(snakemake.input["rgi"]) as tsvfile:
    csvreader = csv.DictReader(tsvfile, delimiter="\t")
    for row in csvreader:
        all_res += extract_rgi(row, aro, [])
            
df = pandas.DataFrame(all_res, columns= ["Software", "Gene", "Resistance Type", "Variant", "Antibiotic resistance prediction"])


writer = pandas.ExcelWriter(snakemake.output["xlsx"])
df.drop_duplicates().to_excel(writer, snakemake.wildcards["sample"], index=False)
writer.save()

df.drop_duplicates().to_csv(snakemake.output["tsv"], index= False, sep="\t")
writer.save()
