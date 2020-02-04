#!/usr/bin/env python

import pandas




def compare_tables(t1,t2, columns):

    RR_count = 0
    SS_count = 0
    RS_count = 0 
    SR_count = 0
    other = 0

    table_a = pandas.read_csv(t1, index_col=0, sep="\t")
    table_b = pandas.read_csv(t2, index_col=0, sep="\t")
    print("file1\tfile2\tsample\tdrug\tsuceptibility1\ttsuceptibility2")
    for sample in table_a.index:
        for column in columns:
            sa = table_a.at[sample, column]
            sb = table_b.at[sample, column]
            if sa != sb:
                print(f"{t1}\t{t2}\t{sample}\t{column}\t{sa}\t{sb}")
            if sa == 'R' and sb == 'R':
                RR_count += 1
            elif sa == 'S' and sb == 'S':
                SS_count += 1
            elif sa == 'R' and sb == 'S':
                RS_count += 1
            elif sa == 'S' and sb == 'R':
                SR_count += 1
            else:
                other+=1


    print(f"RR\t{RR_count}")
    print(f"SS\t{SS_count}")
    print(f"RS\t{RS_count}")
    print(f"SR\t{SR_count}")
    print(f"Other\t{other}")

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-t1", '--table_1', type=str, help="table1")
    parser.add_argument("-t2", '--table_2', type=str, help="table2")

    args = parser.parse_args()

    COLUMNS = ["isoniazid","rifampicin","ethambutol","pyrazinamide"]

    compare_tables(args.table_1, args.table_2, COLUMNS)



