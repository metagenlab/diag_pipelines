#!/usr/bin/env python

def parse_walker_rable(table):
    import re
    with open(table, "r") as f:
        for n, row in enumerate(f):
            data = row.rstrip().split("\t")
            drug = data[0]
            mutation = data[1]	
            details_indel = data[2]
            characterisation = data[3]
            source = data[4]

            
            # skip sensitive mutations
            if characterisation == 'S': 
                continue

            # mutation structure: fabG1_C-15T
            # fabG1	fabG1_T-8C	-8	-8	T	C	SNP
            regex = '(.*)_([A-Za-z]+)([\d\-]+)([A-Za-z]+)'
            s = re.search(regex, mutation)
            
            if s:
                gene = s.group(1)
                REF = s.group(2)
                position = s.group(3)
                ALT = s.group(4)
                if gene == 'rpoB':
                    coli_position = int(position) + 81
                else:
                    coli_position = position
                print(f"{gene}\t{mutation}\t{position}\t{position}\t{REF}\t{ALT}\tSNP")
            else:
                # structure: embA_2723_indel
                # gidB	gidB_141_delC	141	141	delC	delC	INDEL
                regex2 = '(.*)_(\d+)_(.*)'
                s = re.search(regex2, mutation)
                if s:
                    gene = s.group(1)
                    position = s.group(2)
                    change_type = s.group(3)
                    if gene == 'rpoB':
                        coli_position = int(position) + 81
                    else:
                        coli_position = position
                    print(f"{gene}\t{mutation}\t{position}\t{coli_position}\t{change_type}\t{change_type}\tINDEL")
                else:
                    print(mutation)






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