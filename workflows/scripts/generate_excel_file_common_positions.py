import pandas
import itertools
import csv

source_db = {}
source_db["SNP"] = {"CARD": snakemake.input["card"], "Miotto et al. 2017": snakemake.input["miotto"], "Mykrobe": snakemake.input["mykrobe"], "Walker et al. 2015": snakemake.input["walker"]}
nb_sources_snps = len(source_db["SNP"])

#source_db["INDEL"] = {"Mykrobe":"../db/mykrobe_annotated.tsv", "Walker et al. 2015":"../db/walker_resistant_annotated.tsv"}

gene_antibio = pandas.read_csv(snakemake.input['resistance_genes'], sep="\t")
gene_antibio["Antibiotic"] = gene_antibio["Antibiotic"].str.strip()
gene_antibio["Gene"] = gene_antibio["Gene"].str.strip()
genes = sorted(set(list(gene_antibio["Gene"])))

locus_tags = pandas.read_csv(snakemake.input['locus_tag'], sep="\t", index_col=0)

db_combinations = {"4":"four", "3":"three", "2":"two", "1":"one"}

all_databases = {}
for i in source_db["SNP"].keys():
    all_databases[i] = pandas.read_csv(source_db["SNP"][i], sep="\t")
    all_databases[i]["MutatedAminoAcidOrNucleotide"]=all_databases[i]["MutatedAminoAcidOrNucleotide"].str.upper().str.strip()
    all_databases[i]["WildTypeAminoAcidOrNucleotide"]=all_databases[i]["WildTypeAminoAcidOrNucleotide"].str.upper().str.strip()
    all_databases[i]["Position"]=pandas.to_numeric(all_databases[i]["Position"])

common_positions={}
for j in db_combinations.keys():
    common_positions[j] = {}
    for i in set(list(gene_antibio["Antibiotic"])):
        common_positions[j][i] = {}
        common_positions[j][i]["SNP"] = []
        common_positions[j][i]["INDEL"] = []
    
for mutation_type in ["SNP"]:
    all_sources=sorted(list(source_db[mutation_type].keys()))
    for gene in genes:
        pos_each_db={}
        for key in all_sources:
            pos_each_db[key] = set(sorted([int(x) for x in all_databases[key].loc[(all_databases[key]["MutationType"]==mutation_type) & (all_databases[key]["Gene"]==gene), "Position"]]))
        union_4 = sorted(set(pos_each_db[all_sources[0]] & pos_each_db[all_sources[1]] & pos_each_db[all_sources[2]] & pos_each_db[all_sources[3]]))
        common_3 = []
        for db in itertools.combinations(source_db[mutation_type].keys(), 3):
            common_3.append(pos_each_db[db[0]].intersection(pos_each_db[db[1]]).intersection(pos_each_db[db[2]]).difference(union_4))
        union_3 = sorted(set([item for sublist in common_3 for item in sublist]))
        common_2 = []
        for db in itertools.combinations(source_db[mutation_type].keys(), 2):
            common_2.append(pos_each_db[db[0]].intersection(pos_each_db[db[1]]).difference(union_4).difference(union_3))
        union_2 = sorted(set([item for sublist in common_2 for item in sublist]))
        common_1 = []
        for db in source_db[mutation_type].keys():
            common_1.append(pos_each_db[db].difference(union_4).difference(union_3).difference(union_2))
        union_1 = sorted(set([item for sublist in common_1 for item in sublist]))
        unions = {}
        unions['4'] = union_4
        unions['3'] = union_3
        unions['2'] = union_2
        unions['1'] = union_1
        for anti in set(gene_antibio.loc[gene_antibio["Gene"]==gene, "Antibiotic"]):
            for number in ["4", "3", "2", "1"]:
                for pos in unions[number]:
                    for i in source_db[mutation_type].keys():
                        if pos in list(all_databases[i].loc[(all_databases[i]["Gene"]==gene) & (all_databases[i]["MutationType"]==mutation_type), "Position"]):
                            if mutation_type=="SNP":
                                mutated=sorted(set(all_databases[i].loc[(all_databases[i]["MutationType"]==mutation_type) & (all_databases[i]["Gene"]==gene) & (all_databases[i]["Position"]==pos), "MutatedAminoAcidOrNucleotide"]))
                                wildtype=sorted(set(all_databases[i].loc[(all_databases[i]["MutationType"]==mutation_type) & (all_databases[i]["Gene"]==gene) & (all_databases[i]["Position"]==pos), "WildTypeAminoAcidOrNucleotide"]))
                                common_positions[number][anti][mutation_type].append([i, gene, int(pos), ",".join(wildtype), ",".join(mutated)])
                            elif mutation_type=="INDEL":
                                annot=sorted(set(all_databases[i].loc[(all_databases[i]["MutationType"]==mutation_type) & (all_databases[i]["Gene"]==gene) & (all_databases[i]["Position"]==pos), "OriginalAnnotation"]))
                                common_positions[number][anti][mutation_type].append([i, gene, int(pos), annot])

                                
def alternate_gray_background(row, numb):
    if row.empty:
        return row
    ret=row.copy()
    ret[[bool(x) for x in row.index//numb%2]]='background-color:#DCDCDC'
    ret[[not bool(x) for x in row.index//numb%2]]=''
    return(ret)


for j in db_combinations.keys():
    writer = pandas.ExcelWriter(snakemake.output["summary_"+db_combinations[j]])
    all_res = pandas.DataFrame()
    for i in sorted(set(list(gene_antibio["Antibiotic"]))):
        res_panda_snp=pandas.DataFrame(common_positions[j][i]["SNP"], columns=["Source", "Gene", "Position", "WildTypeAminoAcidOrNucleotide", "MutatedAminoAcidOrNucleotide"])
        res_panda_snp.style.apply(alternate_gray_background, axis=None, numb=int(j)).to_excel(writer, i, index=False)

        all_res = all_res.append(res_panda_snp)
    all_res["MutationType"] = "SNP"
    all_res = all_res[~all_res["WildTypeAminoAcidOrNucleotide"].str.contains(",")]
    all_res[["Gene", "Position", "WildTypeAminoAcidOrNucleotide", "MutationType"]].drop_duplicates().to_csv(snakemake.output[db_combinations[j]+"_tsv"], header=True, index=False, sep="\t")
    writer.save()




#    res_panda_indel.style.apply(alternate_gray_background, axis=None).to_excel(writer, "INDEL", index=False)
