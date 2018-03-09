import pandas
import matplotlib
import itertools
from matplotlib_venn import venn3
from matplotlib_venn import venn2

import matplotlib.pyplot as plt

rgi = pandas.read_csv("../db/rgi_annotated_full_2_0_0.csv", sep="\t")
rgi["MutatedAminoAcidOrNucleotide"]=rgi["MutatedAminoAcidOrNucleotide"].str.upper().str.strip()
rgi["WildTypeAminoAcidOrNucleotide"]=rgi["WildTypeAminoAcidOrNucleotide"].str.upper().str.strip()
miotto = pandas.read_csv("../db/miotto_high_moderate_minimum_confidence_annotated.csv", sep="\t")
miotto["MutatedAminoAcidOrNucleotide"]=miotto["MutatedAminoAcidOrNucleotide"].str.upper().str.strip()
miotto["WildTypeAminoAcidOrNucleotide"]=miotto["WildTypeAminoAcidOrNucleotide"].str.upper().str.strip()
mykrobe = pandas.read_csv("../db/mykrobe_annotated.tsv", sep="\t")
mykrobe["MutatedAminoAcidOrNucleotide"]=mykrobe["MutatedAminoAcidOrNucleotide"].str.upper().str.strip()
mykrobe["WildTypeAminoAcidOrNucleotide"]=mykrobe["WildTypeAminoAcidOrNucleotide"].str.upper().str.strip()
walker = pandas.read_csv("../db/walker_resistant_annotated.tsv", sep="\t")
walker["MutatedAminoAcidOrNucleotide"]=walker["MutatedAminoAcidOrNucleotide"].str.upper().str.strip()
walker["WildTypeAminoAcidOrNucleotide"]=walker["WildTypeAminoAcidOrNucleotide"].str.upper().str.strip()


results = [rgi, miotto, mykrobe, walker]
nb_software = len(results)

genes = ["rrs","eis", "rpoB", "rpsL","gyrA", "gyrB", "katG", "pncA", "rrs", "tlyA", "gidB", "inhA", "embA", "embB", "embC", "embR", "ethA", "folC", "iniA", "iniB", "iniC", "kasA", "ndh", "ribD", "thyA", "ahpC", "mabA"]

common_positions = []
for gene in genes:
    for mutation_type in ["SNP","INDEL"]:
        rgi_res = set(sorted([int(x) for x in rgi.loc[(rgi["MutationType"]==mutation_type) & (rgi["Gene"]==gene), "PositionMTB"]]))
        mykrobe_res = set(sorted([int(x) for x in mykrobe.loc[(mykrobe["MutationType"]==mutation_type) & (mykrobe["Gene"]==gene), "PositionMTB"]]))
        walker_res = set(sorted([int(x) for x in walker.loc[(walker["MutationType"]==mutation_type) & (walker["Gene"]==gene), "PositionMTB"]]))
        miotto_res = set(sorted([int(x) for x in miotto.loc[(miotto["MutationType"]==mutation_type) & (miotto["Gene"]==gene), "PositionMTB"]]))
        print(rgi_res)
                         
                            
        common =  sorted(rgi_res & mykrobe_res & miotto_res & walker_res )
        for pos in common:
            mutated_rgi=sorted(set(rgi.loc[(rgi["MutationType"]==mutation_type) & (rgi["Gene"]==gene) & (rgi["PositionMTB"]==pos), "MutatedAminoAcidOrNucleotide"]))
            wildtype_rgi=sorted(set(rgi.loc[(rgi["MutationType"]==mutation_type) & (rgi["Gene"]==gene) & (rgi["PositionMTB"]==pos), "WildTypeAminoAcidOrNucleotide"]))
            common_positions.append(["CARD", gene, str(pos), ",".join(wildtype_rgi), ",".join(mutated_rgi)])

            mutated_mykrobe=sorted(set(mykrobe.loc[(mykrobe["MutationType"]==mutation_type) & (mykrobe["Gene"]==gene) & (mykrobe["PositionMTB"]==pos), "MutatedAminoAcidOrNucleotide"]))
            wildtype_mykrobe=sorted(set(mykrobe.loc[(mykrobe["MutationType"]==mutation_type) & (mykrobe["Gene"]==gene) & (mykrobe["PositionMTB"]==pos), "WildTypeAminoAcidOrNucleotide"]))
            common_positions.append(["Mykrobe", gene, str(pos), ",".join(wildtype_mykrobe), ",".join(mutated_mykrobe)])
            
            mutated_miotto=sorted(set(miotto.loc[(miotto["MutationType"]==mutation_type) & (miotto["Gene"]==gene) & (miotto["PositionMTB"]==pos), "MutatedAminoAcidOrNucleotide"]))
            wildtype_miotto=sorted(set(miotto.loc[(miotto["MutationType"]==mutation_type) & (miotto["Gene"]==gene) & (miotto["PositionMTB"]==pos), "WildTypeAminoAcidOrNucleotide"]))
            common_positions.append(["Miotto", gene, str(pos), ",".join(wildtype_miotto),",".join(mutated_miotto)])
            
            mutated_walker=sorted(set(walker.loc[(walker["MutationType"]==mutation_type) & (walker["Gene"]==gene) & (walker["PositionMTB"]==pos), "MutatedAminoAcidOrNucleotide"]))
            wildtype_walker=sorted(set(walker.loc[(walker["MutationType"]==mutation_type) & (walker["Gene"]==gene) & (walker["PositionMTB"]==pos), "WildTypeAminoAcidOrNucleotide"]))
            common_positions.append(["Walker", gene, str(pos), ",".join(wildtype_walker),",".join(mutated_walker)])

            
def alternate_gray_background(row):
    ret=row.copy()
    ret[[bool(x) for x in row.index//4%2]]='background-color:#DCDCDC'
    ret[[not bool(x) for x in row.index//4%2]]=''
    return(ret)

    
    
res_panda=pandas.DataFrame(common_positions, columns=["Source", "Gene", "Position", "WildTypeAminoAcidorNucleotide", "MutatedAminoAcidOrNucleotide"])

writer = pandas.ExcelWriter("summary.xlsx")
res_panda.style.apply(alternate_gray_background, axis=None).to_excel(writer, index=False)
writer.save()



