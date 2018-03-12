import pandas
import matplotlib
import itertools
from matplotlib_venn import venn3
from matplotlib_venn import venn2

import matplotlib.pyplot as plt

source_db = {}
source_db["SNP"] = {"CARD":"../db/rgi_annotated_full_2_0_0.csv" , "Miotto et al. 2017":"../db/miotto_high_moderate_minimum_confidence_annotated.csv", "Mykrobe":"../db/mykrobe_annotated.tsv", "Walker et al. 2015":"../db/walker_resistant_annotated.tsv"}
source_db["INDEL"] = {"Mykrobe":"../db/mykrobe_annotated.tsv", "Walker et al. 2015":"../db/walker_resistant_annotated.tsv"}

nb_sources = len(source_db)

genes = set(sorted(["rrs","eis", "rpoB", "rpsL","gyrA", "gyrB", "katG", "pncA", "rrs", "tlyA", "gidB", "inhA", "embA", "embB", "embC", "embR", "ethA", "folC", "iniA", "iniB", "iniC", "kasA", "ndh", "ribD", "thyA", "ahpC", "mabA"]))

common_positions = {}
common_positions["SNP"] = []
common_positions["INDEL"] = []

all_res = {}
for i in source_db["SNP"].keys():
    all_res[i] = pandas.read_csv(source_db["SNP"][i], sep="\t")
    all_res[i]["MutatedAminoAcidOrNucleotide"]=all_res[i]["MutatedAminoAcidOrNucleotide"].str.upper().str.strip()
    all_res[i]["WildTypeAminoAcidOrNucleotide"]=all_res[i]["WildTypeAminoAcidOrNucleotide"].str.upper().str.strip()

for mutation_type in ["SNP","INDEL"]:
    for gene in genes:
        key = sorted(source_db[mutation_type].keys())[0]
        common = set(sorted([int(x) for x in all_res[key].loc[(all_res[key]["MutationType"]==mutation_type) & (all_res[key]["Gene"]==gene), "PositionMTB"]]))
        for i in sorted(source_db[mutation_type].keys())[1:]:
            common = common & set(sorted([int(x) for x in all_res[i].loc[(all_res[i]["MutationType"]==mutation_type) & (all_res[i]["Gene"]==gene), "PositionMTB"]]))
        for pos in common:
            for i in source_db[mutation_type].keys():
                if mutation_type=="SNP":
                    mutated=sorted(set(all_res[i].loc[(all_res[i]["MutationType"]==mutation_type) & (all_res[i]["Gene"]==gene) & (all_res[i]["PositionMTB"]==pos), "MutatedAminoAcidOrNucleotide"]))
                    wildtype=sorted(set(all_res[i].loc[(all_res[i]["MutationType"]==mutation_type) & (all_res[i]["Gene"]==gene) & (all_res[i]["PositionMTB"]==pos), "WildTypeAminoAcidOrNucleotide"]))
                    common_positions[mutation_type].append([i, gene, str(pos), ",".join(wildtype), ",".join(mutated)])
                elif mutation_type=="INDEL":
                    annot=sorted(set(all_res[i].loc[(all_res[i]["MutationType"]==mutation_type) & (all_res[i]["Gene"]==gene) & (all_res[i]["PositionMTB"]==pos), "OriginalAnnotation"]))
                    common_positions[mutation_type].append([i, gene, str(pos), annot])
                    

                
def alternate_gray_background(row):
    ret=row.copy()
    ret[[bool(x) for x in row.index//4%2]]='background-color:#DCDCDC'
    ret[[not bool(x) for x in row.index//4%2]]=''
    return(ret)
  
    
res_panda_snp=pandas.DataFrame(common_positions["SNP"], columns=["Source", "Gene", "Position", "WildTypeAminoAcidorNucleotide", "MutatedAminoAcidOrNucleotide"])
res_panda_indel=pandas.DataFrame(common_positions["INDEL"], columns=["Source", "Gene", "Position", "OriginalAnnotation"])

writer = pandas.ExcelWriter("summary.xlsx")
res_panda_snp.style.apply(alternate_gray_background, axis=None).to_excel(writer, "SNP", index=False)
res_panda_indel.style.apply(alternate_gray_background, axis=None).to_excel(writer, "INDEL", index=False)
writer.save()



