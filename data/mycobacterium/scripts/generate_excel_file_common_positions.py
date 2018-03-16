import pandas
import matplotlib
import itertools
from matplotlib_venn import venn3
from matplotlib_venn import venn2

import matplotlib.pyplot as plt

source_db = {}
source_db["SNP"] = {"CARD":"../db/rgi_annotated_full_2_0_0.csv" , "Miotto et al. 2017":"../db/miotto_high_moderate_minimum_confidence_annotated.csv", "Mykrobe":"../db/mykrobe_annotated.tsv", "Walker et al. 2015":"../db/walker_resistant_annotated.tsv"}
source_db["INDEL"] = {"Mykrobe":"../db/mykrobe_annotated.tsv", "Walker et al. 2015":"../db/walker_resistant_annotated.tsv"}

gene_antibio = pandas.read_csv("../db/resistance_genes.csv", sep="\t")
gene_antibio["Antibiotic"]=gene_antibio["Antibiotic"].str.strip()
gene_antibio["Gene"]=gene_antibio["Gene"].str.strip()


nb_sources = len(source_db)

genes = sorted(set(list(gene_antibio["Gene"])))


print(set(list(gene_antibio["Antibiotic"])))

print(genes)

all_res = {}
for i in source_db["SNP"].keys():
    all_res[i] = pandas.read_csv(source_db["SNP"][i], sep="\t")
    all_res[i]["MutatedAminoAcidOrNucleotide"]=all_res[i]["MutatedAminoAcidOrNucleotide"].str.upper().str.strip()
    all_res[i]["WildTypeAminoAcidOrNucleotide"]=all_res[i]["WildTypeAminoAcidOrNucleotide"].str.upper().str.strip()
    all_res[i]["PositionMTB"]=pandas.to_numeric(all_res[i]["PositionMTB"])


common_positions={}
for j in ["4", "3", "2"]:
    common_positions[j] = {}
    for i in set(list(gene_antibio["Antibiotic"])):
        common_positions[j][i] = {}
        common_positions[j][i]["SNP"] = []
        common_positions[j][i]["INDEL"] = []
    
    
for mutation_type in ["SNP"]:
    all_keys=sorted(list(source_db[mutation_type].keys()))
    for gene in genes:
        print(gene)
        pos_each_db={}
        for key in all_keys:
            pos_each_db[key] = set(sorted([int(x) for x in all_res[key].loc[(all_res[key]["MutationType"]==mutation_type) & (all_res[key]["Gene"]==gene), "PositionMTB"]]))
        union_4 = sorted(set(pos_each_db[all_keys[0]] & pos_each_db[all_keys[1]] & pos_each_db[all_keys[2]] & pos_each_db[all_keys[3]]))
        common_3 = []
        for db in itertools.combinations(source_db[mutation_type].keys(), 3):
            common_3.append(pos_each_db[db[0]].intersection(pos_each_db[db[1]]).intersection(pos_each_db[db[2]]).difference(union_4))
        union_3 = sorted(set([item for sublist in common_3 for item in sublist]))
        common_2 = []
        for db in itertools.combinations(source_db[mutation_type].keys(), 2):
            common_2.append(pos_each_db[db[0]].intersection(pos_each_db[db[1]]).difference(union_4).difference(union_3))
        union_2 = sorted(set([item for sublist in common_2 for item in sublist]))
        unions = {}
        unions['4'] = union_4
        unions['3'] = union_3
        unions['2'] = union_2
        for anti in set(gene_antibio.loc[gene_antibio["Gene"]==gene, "Antibiotic"]):
            for number in ["4", "3", "2"]:
                for pos in unions[number]:
                    for i in source_db[mutation_type].keys():
                        if pos in list(all_res[i].loc[(all_res[i]["Gene"]==gene) & (all_res[i]["MutationType"]==mutation_type), "PositionMTB"]):
                            if mutation_type=="SNP":
                                mutated=sorted(set(all_res[i].loc[(all_res[i]["MutationType"]==mutation_type) & (all_res[i]["Gene"]==gene) & (all_res[i]["PositionMTB"]==pos), "MutatedAminoAcidOrNucleotide"]))
                                wildtype=sorted(set(all_res[i].loc[(all_res[i]["MutationType"]==mutation_type) & (all_res[i]["Gene"]==gene) & (all_res[i]["PositionMTB"]==pos), "WildTypeAminoAcidOrNucleotide"]))
                                common_positions[number][anti][mutation_type].append([i, gene, str(pos), ",".join(wildtype), ",".join(mutated)])
                            elif mutation_type=="INDEL":
                                annot=sorted(set(all_res[i].loc[(all_res[i]["MutationType"]==mutation_type) & (all_res[i]["Gene"]==gene) & (all_res[i]["PositionMTB"]==pos), "OriginalAnnotation"]))
                                common_positions[number][anti][mutation_type].append([i, gene, str(pos), annot])

def alternate_gray_background(row, numb):
    if row.empty:
        return row
    ret=row.copy()
    ret[[bool(x) for x in row.index//numb%2]]='background-color:#DCDCDC'
    ret[[not bool(x) for x in row.index//numb%2]]=''
    return(ret)
  
for j in ["4", "3", "2"]:
    writer = pandas.ExcelWriter("summary_common_to_"+j+"_dbs.xlsx")
    for i in sorted(set(list(gene_antibio["Antibiotic"]))):
        res_panda_snp=pandas.DataFrame(common_positions[j][i]["SNP"], columns=["Source", "Gene", "Position", "WildTypeAminoAcidorNucleotide", "MutatedAminoAcidOrNucleotide"])
        res_panda_snp.style.apply(alternate_gray_background, axis=None, numb=int(j)).to_excel(writer, i, index=False)
    writer.save()
#    res_panda_indel.style.apply(alternate_gray_background, axis=None).to_excel(writer, "INDEL", index=False)
