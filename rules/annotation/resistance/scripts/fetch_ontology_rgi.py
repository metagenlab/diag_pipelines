import pandas
import pronto
import csv
import re

aro = pronto.Ontology(snakemake.params["ontology_aro"])
ro = pronto.Ontology(snakemake.params["ontology_ro"])
mo = pronto.Ontology(snakemake.params["ontology_mo"])


aro.merge(ro)
aro.merge(mo)

ontology_results = []

rgi_results = pandas.read_csv(snakemake.input['rgi'], sep="\t")

cat = snakemake.params["aro_categories"]
categories = pandas.read_csv(cat, sep="\t")

resistance_mechanism = categories.loc[categories['ARO Category'] == 'Resistance Mechanism']
obo2resistance_mechanism = pandas.Series(resistance_mechanism['ARO Name'].values, index=resistance_mechanism['ARO Accession']).to_dict()
AMR_family = categories.loc[categories['ARO Category'] == 'AMR Gene Family']
obo2AMR_family = pandas.Series(AMR_family['ARO Name'].values, index=AMR_family['ARO Accession']).to_dict()
drug_class = categories.loc[categories['ARO Category'] == 'Drug Class']
obo2drug_class = pandas.Series(drug_class['ARO Name'].values, index=drug_class['ARO Accession']).to_dict()

def get_drug_class(obo_drug, obo2drug_class):
    if obo_drug.id in obo2drug_class:
        return obo2drug_class[obo_drug.id]
    for parent in obo_drug.rparents():
        if parent.id in obo2drug_class:
            return obo2drug_class[parent.id]


def get_AMR_family(obo_gene, obo2AMR_family):
    if obo_gene.id in obo2AMR_family:
        return obo2AMR_family[obo_gene.id]
    for parent in obo_gene.rparents():
        if parent.id in obo2AMR_family:
            return obo2AMR_family[parent.id]

def get_all_path(obo_list, all_obo):
    new_obo = []
    for obo in obo_list:
        new_obo+=obo.parents
        all_obo+=obo.parents
    if len(new_obo) > 0:
        return get_all_path(new_obo, all_obo=all_obo)
    else:
        return all_obo

def search_participating_mechanism(obo_term, mechanism_list):
    for relation in obo_term.relations:
        if relation.obo_name=="participates_in":
            if obo_term.relations[relation][0].id in mechanism_list:
                return obo_term.relations[relation][0]
            else:
                # search mechanism in antibiotic path
                all_path = get_all_path(obo_term.relations[relation], [])
                try:
                    return [i for i in all_path if i.id in mechanism_list][0]
                except IndexError:
                    continue
        elif relation.obo_name in ["regulates"]:
            # classify as regulator
            relation.name = 'regulator'
            return (relation)

for index, row in rgi_results.iterrows():
    if row["Model_type"] == "protein variant model":
        gene = row["Best_Hit_ARO"]
        aro_term = "ARO:%s" % row["ARO"]
        snp = row["SNPs_in_Best_Hit_ARO"]

        term = aro[aro_term]
        fam = get_AMR_family(term, obo2AMR_family)
        for parent in term.rparents():
            mechanism = search_participating_mechanism(parent, obo2resistance_mechanism.keys())
            if mechanism:
                break
        print (snp)
        pos = re.sub("[^0-9]", "", snp)
        ref = snp[0]
        mut = snp[-1]
        for child in term.relations:
            if child.obo_name=="confers_resistance_to_drug" or child.obo_name=="confers_resistance_to":
                for antibiotic in term.relations[child]:
                    anti = antibiotic.name.replace("antibiotic", "").strip()
                    anti_class = get_drug_class(antibiotic, obo2drug_class)
                    if anti_class is not None:
                        anti_class = anti_class.replace("antibiotic", "")
                    ontology_results.append([term.id, term.name, "rgi", gene, "variant model", ref+str(pos)+mut, anti, anti_class, mechanism.name, fam])

    elif row["Model_type"] == "protein homolog model":
        hit_list = []

        term = aro["ARO:%s" % row["ARO"]]
        print(term, term.id, term.name, '---------------')
        fam = get_AMR_family(term, obo2AMR_family)

        # search recursively
        for parent in term.rparents():
            mechanism = search_participating_mechanism(parent, obo2resistance_mechanism.keys())
            if mechanism:
                break
        gene = row["Best_Hit_ARO"].replace(" ", "_").strip()
        print('relations', term.relations)
        for child in term.relations:
            print(child.obo_name)
            if child.obo_name=="confers_resistance_to_drug" or child.obo_name=="confers_resistance_to":
                for antibiotic in term.relations[child]:
                    anti_class = get_drug_class(antibiotic, obo2drug_class)
                    if anti_class is not None:
                        anti_class = anti_class.replace("antibiotic", "")
                    anti = antibiotic.name.replace("antibiotic", "").strip()

                    # check if it is a parent of already included antibio
                    # if yes, skip (same thing)
                    parent_of = False
                    for hit in hit_list:
                        if child in hit.rparents():
                            parent_of = True
                            break
                    if parent_of:
                        continue
                    if mechanism is not None:
                        mechanism_name = mechanism.name
                    else:
                        mechanism_name = 'Unknown'

                    print('appaned:', [term.id, term.name, "rgi", gene, "homology model", "", anti, anti_class, mechanism_name, fam])
                    ontology_results.append([term.id, term.name, "rgi", gene, "homology model", "", anti, anti_class, mechanism_name, fam])
            print('parents', term.rparents())
            for parent in term.rparents():
                for child in parent.relations:
                    print(child.obo_name)
                    if child.obo_name=="confers_resistance_to_drug" or child.obo_name=="confers_resistance_to":
                        for antibiotic in parent.relations[child]:
                            anti_class = get_drug_class(antibiotic, obo2drug_class)
                            anti = antibiotic.name.replace("antibiotic", "").strip()

                            # check if it is a parent of already included antibio
                            # if yes, skip (same thing)
                            parent_of = False
                            for hit in hit_list:
                                if parent in hit.rparents():
                                    parent_of=True
                                    break
                            if parent_of:
                                # sub antibiotic class already into list, skip
                                continue
                            if anti_class is not None:
                                anti_class = anti_class.replace("antibiotic", "")
                            if mechanism is not None:
                                mechanism_name = mechanism.name
                            else:
                                mechanism_name = 'Unknown'
                            print('append:', [term.id, term.name, "rgi", gene, "homology model", "", anti, anti_class, mechanism_name, fam])
                            ontology_results.append([term.id, term.name, "rgi", gene, "homology model", "", anti, anti_class, mechanism_name, fam])
                            hit_list.append(parent)


df = pandas.DataFrame(ontology_results, columns= ["ARO","Name", "Software", "Gene", "Resistance Type", "Variant", "Antibiotic resistance prediction", "Class", "Mechanism", "AMR family"])

writer = pandas.ExcelWriter(snakemake.output["xlsx"])
df.drop_duplicates().to_excel(writer, snakemake.wildcards["sample"], index=False)
writer.save()

df.drop_duplicates().to_csv(snakemake.output["tsv"], index= False, sep="\t")
writer.save()
