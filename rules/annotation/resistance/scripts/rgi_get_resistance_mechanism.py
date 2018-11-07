import pandas
import pronto
import csv
import re
from pronto import Relationship

aro = pronto.Ontology(snakemake.params["ontology_aro"])
ro = pronto.Ontology(snakemake.params["ontology_ro"])
mo = pronto.Ontology(snakemake.params["ontology_mo"])
cat = snakemake.params["aro_categories"]
categories = pandas.read_csv(cat, sep="\t")

resistance_mechanism = categories.loc[categories['ARO Category'] == 'Resistance Mechanism']

obo2resistance_mechanism = pandas.Series(resistance_mechanism['ARO Name'].values, index=resistance_mechanism['ARO Accession']).to_dict()

aro.merge(ro)
aro.merge(mo)

ontology_results = []

rgi_results = pandas.read_csv(snakemake.input['rgi'], sep="\t")


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
                all_path = get_all_path(obo_term.relations[relation], [])
                try:
                    return [i for i in all_path if i.id in mechanism_list][0]
                except IndexError:
                    continue
        elif relation.obo_name in ["regulates"]:
            # classify as regulator
            relation.name = 'regulator'
            return (relation)

term_nr_list = []
for index, row in rgi_results.iterrows():
    gene = row["Best_Hit_ARO"]
    aro_id = str(row["ARO"])
    aro_term = "ARO:%s" % row["ARO"]
    term = aro[aro_term]
    for i in term.rparents():
        mechanism = search_participating_mechanism(i, obo2resistance_mechanism.keys())
        if mechanism:
            break
    # remove redundancy
    if mechanism is not None:
        if [gene, mechanism.name] not in term_nr_list:
            term_nr_list.append([gene, mechanism.name])
    else:
        term_nr_list.append([gene, "Unknown"])

df = pandas.DataFrame(term_nr_list, columns= ["Gene", "Resistance Mechanism"])


writer = pandas.ExcelWriter(snakemake.output["xlsx"])
df.drop_duplicates().to_excel(writer, snakemake.wildcards["sample"], index=False)
writer.save()

df.drop_duplicates().to_csv(snakemake.output["tsv"], index= False, sep="\t")
writer.save()
