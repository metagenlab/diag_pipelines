#!/usr/bin/env python


def fasta_list2bins(fasta_list):
    from Bio import SeqIO
    import sys
    import os
    import re
    
    for fasta in fasta_list:
        with open(fasta, 'r') as f:
            records = SeqIO.parse(f, 'fasta')
            for record in records:
                sys.stdout.write("%s\t%s\n" % (record.name,
                                               re.sub('\.','_',os.path.basename(fasta.split(".fasta")[0]))))
    
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_fasta', type=str, help="input fasta files", nargs="+")

    args = parser.parse_args()
    fasta_list2bins(args.input_fasta)
    
