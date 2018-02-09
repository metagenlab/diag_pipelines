import pandas
import pronto
import csv
import re



out_csv = snakemake.output[0]
out_xlsx = snakemake.output[1]

aro = pronto.Ontology(snakemake.params["ontology_aro"])
ro = pronto.Ontology(snakemake.params["ontology_ro"])
mo = pronto.Ontology(snakemake.params["ontology_mo"])

currated_genes = snakemake.params["currated_genes"]

aro.merge(ro)
aro.merge(mo)

all_res = []


def extract_rgi(row, aro_ont, gene_list):
    res = []
    if row["SNP"] != "n/a":
        best_hit = row["Best_Hit_ARO"].split()
        genes = list(set(best_hit))
        if len(genes) == 1:
            gene = genes[0]
            terms = row["ARO"].split(",")
            SNPS = row["SNP"].split(",")
            for i in range(len(terms)):
                term = aro_ont[terms[i].strip()]
                snp = SNPS[i].strip()
                pos = re.sub("[^0-9]", "", snp)
                ref = snp[0]
                mut = snp[-1]
                for child in term.relations:
                    if child.obo_name=="confers_resistance_to_drug" or  child.obo_name=="confers_resistance_to":
                        for antibiotic in term.relations[child]:
                            anti = antibiotic.name.replace("antibiotic", "").strip()
                            res.append(["rgi", gene, "variant", ref+str(pos)+mut, anti])
        else:
            raise Exception("Problem parsing the predicted ARO Model from RGI")
    elif "intrinsic" not in row["Best_Hit_ARO"]:
        term = aro_ont[row["ARO"]]
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


def extract_mykrobe(row, aro_ont=None, gene_list=None):
    res = []
    if row["susceptibility"] == "R" and row["variants (gene:alt_depth:wt_depth:conf)"] =="":
        for result in  row["genes (prot_mut-ref_mut:percent_covg:depth)"].split(";"):
            gene = result.split(":")[0]
            anti = row["drug"].lower()
            res.append(["mykrobe", gene, "gene presence", "", anti])
    elif row["susceptibility"] == "R" and row["variants (gene:alt_depth:wt_depth:conf)"] !="" and "del" not in row["variants (gene:alt_depth:wt_depth:conf)"] and "ins" not in row["variants (gene:alt_depth:wt_depth:conf)"]:
        anti = row["drug"].lower()
        for result in  row["variants (gene:alt_depth:wt_depth:conf)"].split(";"):
            gene = result.split("_")[0]
            mutation = result.split('_')[1].split("-")[0]
            ref = mutation[0]
            mut = mutation[-1]
            pos = re.sub("[^0-9]", "", mutation)
            res.append(["mykrobe", gene, "variant", ref+str(pos)+mut, anti])
    elif "del" in row["variants (gene:alt_depth:wt_depth:conf)"] or "ins" in row["variants (gene:alt_depth:wt_depth:conf)"]:
        raise Exception("Insertion or deletion confering resistance have been detected by Mykrobe in {0}, parsing currently not implemented".format(sample))
    return(res)


def extract_abricate(row, aro_ont=None, gene_list=None):
    res = []
    gene = row["GENE"]
    res.append(["abricate", gene, "gene presence", "", ""])
    return(res)

for file_tsv in snakemake.input:
    if 'rgi' in file_tsv:
        func = extract_rgi
    elif 'mykrobe' in file_tsv:
        func = extract_mykrobe
    elif 'abricate' in file_tsv:
        func = extract_abricate
    with open(file_tsv, "r") as tsvfile:
        csvreader = csv.DictReader(tsvfile, delimiter="\t")
        for row in csvreader:
            all_res += func(row, aro, currated_genes)
        


            
df = pandas.DataFrame(all_res, columns= ["Software", "Gene", "Resistance Type", "Variant", "Antibiotic resistance prediction"])

writer = pandas.ExcelWriter(out_xlsx)
df.to_excel(writer, snakemake.wildcards["sample"], index=False)
writer.save()

df.to_csv(out_csv, index= False, sep="\t")
writer.save()
