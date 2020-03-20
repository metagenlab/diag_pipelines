#!/usr/bin/env python


GENE_LIST = ["tlyA", 
             "thyA", 
             "rpsL",
             "rpsA",
             "rpoB",
             "ribD",
             "pncA",
             "ndh",
             "katG",
             "kasA",
             "iniC",
             "iniB",
             "iniA",
             "inhA",
             "gyrB",
             "gyrA",
             "gidB",
             "folC",
             "embR",
             "embB",
             "embA",
             "ethA"]


def parse_card_json(json_file):
    import re
    import json

    target_models = ["protein variant model", "rRNA gene variant model", ""]
    

    with open(json_file, "r") as read_file:
        card_data = json.load(read_file)

        for model in card_data:
            antibiotic = []
            drug_class = None
            if "_" in model:
                continue
            if card_data[model]['model_type'] == "protein variant model":
                #print(card_data[model]["model_sequences"]["sequence"])
                sequence_list = card_data[model]["model_sequences"]["sequence"]
                for sequence in sequence_list:
                    if "Mycobacterium tuberculosis" in card_data[model]["model_sequences"]["sequence"][sequence]["NCBI_taxonomy"]["NCBI_taxonomy_name"]:

                        snp_data = card_data[model]['model_param']['snp']

                        for param in card_data[model]['model_param']:

                            #if card_data[model]['model_param'][param]['param_type'] == 'multiple resistance variants':
                            #    print("multiple resistance variants",card_data[model]['model_param'][param])
                            #    print()

                            if param != 'blastp_bit_score':
                                param_type = card_data[model]['model_param'][param]['param_type']
                                param_values = card_data[model]['model_param'][param]["param_value"]
                                model_name = card_data[model]["model_name"].split(" ")[2]

                                if "mutations conferring resistance to" in card_data[model]["model_name"]:
                                    s = re.search("(.*) (.*) mutations conferring resistance to (.*)", card_data[model]["model_name"])
                                elif "with mutation conferring resistance to" in card_data[model]["model_name"]:
                                    s = re.search("(.*) (.*) with mutation conferring resistance to (.*)", card_data[model]["model_name"])
                                elif "mutant conferring resistance to" in card_data[model]["model_name"]:
                                    s = re.search("(.*) (.*) mutant conferring resistance to (.*)", card_data[model]["model_name"])
                                elif "mutants conferring resistance to" in card_data[model]["model_name"]:
                                    s = re.search("(.*) (.*) mutants conferring resistance to (.*)", card_data[model]["model_name"])
                                elif "mutation conferring resistance to" in card_data[model]["model_name"]:
                                    s = re.search("(.*) (.*) mutation conferring resistance to (.*)", card_data[model]["model_name"])                                    
                                elif "conferring resistance to" in card_data[model]["model_name"]:
                                    s = re.search("(.*) (.*) conferring resistance to (.*)", card_data[model]["model_name"])
                                else:
                                    raise IOError("unexpeted model name: %s" % card_data[model]["model_name"])

                                gene = s.group(2)

                                
                                aro = card_data[model]["ARO_accession"]

                                for i in card_data[model]["ARO_category"]:
                                    if (card_data[model]["ARO_category"][i]["category_aro_class_name"] == "Drug Class"):
                                        drug_class = card_data[model]["ARO_category"][i]["category_aro_name"]
                                    if (card_data[model]["ARO_category"][i]["category_aro_class_name"] == "Antibiotic"):
                                        antibiotic.append(card_data[model]["ARO_category"][i]["category_aro_name"])

                                if len(antibiotic) == 0:
                                    category = drug_class
                                # if resistance to multiple antibiotics of the same class, report the class
                                elif len(set(antibiotic)) > 1:
                                    category = re.sub(" antibiotic","", drug_class)
                                # if resistant to a single antibiotic, report the antibiotic name
                                else:
                                    category = re.sub(" antibiotic","", antibiotic[0])
                                for key in param_values:
                                    value = card_data[model]['model_param'][param]["param_value"][key]
                                    
                                    print(f"{aro}\t{gene}\t{param_type}\t{value}\t{category}")
                                




'''

{'param_type': 'single resistance variant', 'param_description': 'A nucleotide or amino acid substitution that confers elevated resistance to antibiotic(s) relative to wild type. The most common type encoded in the CARD is an amino acid substitution gleaned from the literature with format [wild-type][position][mutation], e.g. R184Q. When present in the associated gene or protein, a single resistance variant confers resistance to an antibiotic drug or drug class. Single resistance variants are used by the protein variant and rRNA mutation models to detect antibiotic resistance from submitted sequences.', 'param_type_id': '36301', 'param_value': {'8669': 'D426N', '8670': 'R447L', '8671': 'E466V', '8672': 'S416A', '8673': 'D426V', '8674': 'R447K', '8675': 'E466K', '8814': 'R377G', '8815': 'R389P', '8816': 'E399K', '8817': 'D409N', '8818': 'V423F', '8819': 'R457T', '8820': 'D465Y', '8834': 'S366A', '8835': 'S366V', '8836': 'S464T', '8837': 'I139R', '8838': 'V130I', '8846': 'K444F'}, 'clinical': {'8669': 'D426N', '8670': 'R447L', '8671': 'E466V', '8672': 'S416A', '8673': 'D426V', '8674': 'R447K', '8675': 'E466K', '8814': 'R377G', '8815': 'R389P', '8816': 'E399K', '8817': 'D409N', '8818': 'V423F', '8819': 'R457T', '8820': 'D465Y', '8834': 'S366A', '8835': 'S366V', '8836': 'S464T', '8837': 'I139R', '8838': 'V130I', '8846': 'K444F'}}

model_id
model_name
model_type
model_type_id
model_description
model_param
model_sequences
ARO_accession
ARO_id
ARO_name
ARO_description
ARO_category
'''



header = ["Gene",
          "OriginalAnnotation",
          "Position",
          "PositionAlternative",
          "WildTypeAminoAcidOrNucleotide",
          "MutatedAminoAcidOrNucleotide",
          "MutationType"]


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", '--file', type=str, help="input CARD json file")


    args = parser.parse_args()

    parse_card_json(args.file)