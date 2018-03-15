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

common_positions={}
for i in set(list(gene_antibio["Antibiotic"])):
    common_positions[i] = {}
    common_positions[i]["SNP"] = []
    common_positions[i]["INDEL"] = []

all_res = {}
for i in source_db["SNP"].keys():
    all_res[i] = pandas.read_csv(source_db["SNP"][i], sep="\t")
    all_res[i]["MutatedAminoAcidOrNucleotide"]=all_res[i]["MutatedAminoAcidOrNucleotide"].str.upper().str.strip()
    all_res[i]["WildTypeAminoAcidOrNucleotide"]=all_res[i]["WildTypeAminoAcidOrNucleotide"].str.upper().str.strip()




common_positions_3={}

for i in set(list(gene_antibio["Antibiotic"])):
    common_positions_3[i] = {}
    common_positions_3[i]["SNP"] = []
    common_positions_3[i]["INDEL"] = []
    
for mutation_type in ["SNP"]:
    all_keys=sorted(list(source_db[mutation_type].keys()))
    for gene in genes:
        print(gene)
        pos_each_db={}
        for key in all_keys:
            pos_each_db[key] = set(sorted([int(x) for x in all_res[key].loc[(all_res[key]["MutationType"]==mutation_type) & (all_res[key]["Gene"]==gene), "PositionMTB"]]))
        common_4 = pos_each_db[all_keys[0]] & pos_each_db[all_keys[1]] & pos_each_db[all_keys[2]] & pos_each_db[all_keys[3]] 
        for anti in set(gene_antibio.loc[gene_antibio["Gene"]==gene, "Antibiotic"]):
            for pos in sorted(common_4):
                for i in source_db[mutation_type].keys():
                    if mutation_type=="SNP":
                        mutated=sorted(set(all_res[i].loc[(all_res[i]["MutationType"]==mutation_type) & (all_res[i]["Gene"]==gene) & (all_res[i]["PositionMTB"]==pos), "MutatedAminoAcidOrNucleotide"]))
                        wildtype=sorted(set(all_res[i].loc[(all_res[i]["MutationType"]==mutation_type) & (all_res[i]["Gene"]==gene) & (all_res[i]["PositionMTB"]==pos), "WildTypeAminoAcidOrNucleotide"]))
                        common_positions[anti][mutation_type].append([i, gene, str(pos), ",".join(wildtype), ",".join(mutated)])
                    elif mutation_type=="INDEL":
                        annot=sorted(set(all_res[i].loc[(all_res[i]["MutationType"]==mutation_type) & (all_res[i]["Gene"]==gene) & (all_res[i]["PositionMTB"]==pos), "OriginalAnnotation"]))
                        common_positions[anti][mutation_type].append([i, gene, str(pos), annot])
        common_3 = []
        for db in itertools.combinations(source_db[mutation_type].keys(), 3):
            common_3.append(pos_each_db[db[0]] & pos_each_db[db[1]] & pos_each_db[db[2]])
        union_3 = common_3[0].union(common_3[1].union(common_3[2])).difference(common_4)
        for pos in sorted(union_3):
            for i in source_db[mutation_type].keys():
#                print(list(all_res[i].loc[all_res[i]["Gene"]==gene, "PositionMTB"]))
                if pos in list(all_res[i].loc[all_res[i]["Gene"]==gene, "PositionMTB"]):
                    if mutation_type=="SNP":
                        mutated=sorted(set(all_res[i].loc[(all_res[i]["MutationType"]==mutation_type) & (all_res[i]["Gene"]==gene) & (all_res[i]["PositionMTB"]==pos), "MutatedAminoAcidOrNucleotide"]))
                        wildtype=sorted(set(all_res[i].loc[(all_res[i]["MutationType"]==mutation_type) & (all_res[i]["Gene"]==gene) & (all_res[i]["PositionMTB"]==pos), "WildTypeAminoAcidOrNucleotide"]))
                        common_positions_3[anti][mutation_type].append([i, gene, str(pos), ",".join(wildtype), ",".join(mutated)])
                    elif mutation_type=="INDEL":
                        annot=sorted(set(all_res[i].loc[(all_res[i]["MutationType"]==mutation_type) & (all_res[i]["Gene"]==gene) & (all_res[i]["PositionMTB"]==pos), "OriginalAnnotation"]))
                        common_positions_3[anti][mutation_type].append([i, gene, str(pos), annot])
                    

        
            
                
def alternate_gray_background3(row):
    if row.empty:
        return row
    ret=row.copy()
    ret[[bool(x) for x in row.index//3%2]]='background-color:#DCDCDC'
    ret[[not bool(x) for x in row.index//3%2]]=''
    return(ret)
  
                    

                
def alternate_gray_background(row):
    if row.empty:
        return row
    ret=row.copy()
    ret[[bool(x) for x in row.index//4%2]]='background-color:#DCDCDC'
    ret[[not bool(x) for x in row.index//4%2]]=''
    return(ret)
  
writer = pandas.ExcelWriter("summary.xlsx")

for i in sorted(set(list(gene_antibio["Antibiotic"]))):
    print(i)
    res_panda_snp=pandas.DataFrame(common_positions[i]["SNP"], columns=["Source", "Gene", "Position", "WildTypeAminoAcidorNucleotide", "MutatedAminoAcidOrNucleotide"])

    #    res_panda_indel=pandas.DataFrame(common_positions["INDEL"], columns=["Source", "Gene", "Position", "OriginalAnnotation"])



    res_panda_snp.style.apply(alternate_gray_background, axis=None).to_excel(writer, i, index=False)
writer.save()
#    res_panda_indel.style.apply(alternate_gray_background, axis=None).to_excel(writer, "INDEL", index=False)




  
writer = pandas.ExcelWriter("summary3.xlsx")

for i in sorted(set(list(gene_antibio["Antibiotic"]))):
    print(i)
    res_panda_snp=pandas.DataFrame(common_positions_3[i]["SNP"], columns=["Source", "Gene", "Position", "WildTypeAminoAcidorNucleotide", "MutatedAminoAcidOrNucleotide"])

    #    res_panda_indel=pandas.DataFrame(common_positions["INDEL"], columns=["Source", "Gene", "Position", "OriginalAnnotation"])



    res_panda_snp.style.apply(alternate_gray_background3, axis=None).to_excel(writer, i, index=False)
writer.save()
#    res_panda_indel.style.apply(alternate_gray_background, axis=None).to_excel(writer, "INDEL", index=False)



