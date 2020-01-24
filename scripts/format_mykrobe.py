#!/usr/bin/env python

def parse_walker_rable(json_file):
    import re
    import json
    with open(json_file, "r") as read_file:
        data = json.load(read_file)
        for entry in data:
        
            # mutation structure: rpoB_TGA1289T
            # rpoB_TGA1328T
            # rpoB_TGACCCACAA1328GGCCCCA
            
            regex = '(.*)_([A-Za-z]+)([\d\-]+)([A-Za-z]+)'
            s = re.search(regex, entry)
            if s:
                gene = s.group(1)
                REF = s.group(2)
                position = s.group(3)
                ALT = s.group(4)
                if gene == 'rpoB':
                    coli_position = int(position) + 81
                else:
                    coli_position = position
                if (len(REF) > 1) or (len(ALT) > 1):
                    change_type = 'INDEL'
                else:
                    change_type = 'SNP'
                print(f"{gene}\t{entry}\t{position}\t{coli_position}\t{REF}\t{ALT}\t{change_type}")
            else:
                print("PROBLEM", entry)






header = ["Gene",
          "OriginalAnnotation",
          "Position",
          "PositionAlternative",
          "WildTypeAminoAcidOrNucleotide",
          "MutatedAminoAcidOrNucleotide",
          "MutationType"]


### 13 genes ###
# ahpC
# eis
# embB
# fabG1
# Gene
# gidB
# gyrA
# inhA
# atG
# pncA
# rpoB
# rpsA
# rpsL
# rrs
# tlyA

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", '--file', type=str, help="input walker tsv file")


    args = parser.parse_args()

    parse_walker_rable(args.file)