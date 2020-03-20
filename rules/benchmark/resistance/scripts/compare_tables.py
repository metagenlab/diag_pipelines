#!/usr/bin/env python

import pandas

def compare_tables(t1,
                   t2, 
                   columns, 
                   disrepancies,
                   FP,
                   FN,
                   counts,
                   stats):

    # known synonyms for CARD
    # RGI: rifampicin = rifamycin
    # RGI: ethambutol = polyamine  
    SYNONYMS = {
        "rifampicin" : "rifamycin",
        "ethambutol" : "polyamine"
    }

    table_a = pandas.read_csv(t1, index_col=0, sep="\t")
    table_b = pandas.read_csv(t2, index_col=0, sep="\t")


    disrepancies_f = open(disrepancies, "w")
    FP_f = open(FP, "w")
    FN_f = open(FN, "w")
    counts_f = open(counts, "w")
    stats_f = open(stats, "w")

    ref = t1.split(".")[0]
    tool = t2.split("/")[-1].split("_benchmark.tsv")[0]

    disrepancies_f.write("file1\tfile2\tsample\tdrug\tsuceptibility1\tsuceptibility2\n")
    antibiotic2counts = {} 
    for antibiotic in columns:
        antibiotic2counts[antibiotic] = {"reference" : { "R" : 0, "S": 0, "other": 0}  , "tool" :  { "R" : 0, "S": 0, "other": 0}}  
        antibiotic2counts[antibiotic]["RR"] = 0
        antibiotic2counts[antibiotic]["SS"] = 0
        antibiotic2counts[antibiotic]["RS"] = 0
        antibiotic2counts[antibiotic]["SR"] = 0
        antibiotic2counts[antibiotic]["other"] = 0

    for sample in table_a.index:
        for antibiotic in columns:

            sa = table_a.at[sample, antibiotic]
            colref = antibiotic
            if antibiotic not in table_b.columns.tolist():
                try:
                    antibiotic = SYNONYMS[antibiotic]
                except KeyError:
                    raise IOError("Unkown antibiotic:", antibiotic)

            sb = table_b.at[sample, antibiotic]

            if sa == 'R':
                antibiotic2counts[colref]["reference"]["R"] +=1 
            elif sa == 'S':
                antibiotic2counts[colref]["reference"]["S"] +=1
            elif sa not in ["S", "R"] :
                antibiotic2counts[colref]["reference"]["other"] +=1
            if sb == 'R':
                antibiotic2counts[colref]["tool"]["R"] +=1 
            elif sb == 'S':
                antibiotic2counts[colref]["tool"]["S"] +=1
            elif sb not in ["S", "R"] :
                antibiotic2counts[colref]["tool"]["other"] +=1


            if sa != sb:
                disrepancies_f.write(f"{ref}\t{tool}\t{sample}\t{antibiotic}\t{sa}\t{sb}\n")
            if sa == 'R' and sb == 'R':
                antibiotic2counts[colref]["RR"] += 1
            elif sa == 'S' and sb == 'S':
                antibiotic2counts[colref]["SS"] += 1
            elif sa == 'R' and sb == 'S':
                antibiotic2counts[colref]["RS"] += 1
                # false negative
                FN_f.write(f"{tool}\t{sample}\t{antibiotic}\n")
            elif sa == 'S' and sb == 'R':
                antibiotic2counts[colref]["SR"] += 1
                # false positive
                FP_f.write(f"{tool}\t{sample}\t{antibiotic}\n")
            else:
                antibiotic2counts[colref]["other"] += 1
    

    counts_f.write("antibiotic\tR_ref\tS_ref\tOther_ref\tR_tool\tS_tool\tOther_tool\n")
    for antibiotic in antibiotic2counts:
        # detail
        counts_f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (antibiotic,
                                                         antibiotic2counts[antibiotic]["reference"]["R"],
                                                         antibiotic2counts[antibiotic]["reference"]["S"],
                                                         antibiotic2counts[antibiotic]["reference"]["other"],
                                                         antibiotic2counts[antibiotic]["tool"]["R"],
                                                         antibiotic2counts[antibiotic]["tool"]["S"],
                                                         antibiotic2counts[antibiotic]["tool"]["other"]))

    # comparison
    stats_f.write("\nantibiotic\tRR\tSS\tRS\tSR\tother\tTOTAL\tsensitivity\tspecificity\tprecision\n")
    for antibiotic in antibiotic2counts:

        TP = antibiotic2counts[antibiotic]["RR"]
        TN = antibiotic2counts[antibiotic]["SS"]
        FP = antibiotic2counts[antibiotic]["SR"]
        FN = antibiotic2counts[antibiotic]["RS"]

        # sensitivity (recall)
        # TP / (TP+FN)
        # fraction of resistant strains predicted as resistant

        sensitivity = round(TP / (TP + FN), 3)

        # specificity
        # TN / (TN + FP)
        # fraction fo suceptible strain prediected as suceptible

        specificity = round(TN / (TN + FP), 3)

        # precision
        # TP / (TP + FP)
        # fraction of resistant strains predicted as resistant

        precision = round(TP / (TP + FP), 3)

        s = antibiotic2counts[antibiotic]["RR"] 
        s += antibiotic2counts[antibiotic]["SS"]
        s += antibiotic2counts[antibiotic]["RS"]
        s += antibiotic2counts[antibiotic]["SR"]
        s += antibiotic2counts[antibiotic]["other"]
        stats_f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (antibiotic,
                                                                    antibiotic2counts[antibiotic]["RR"],
                                                                    antibiotic2counts[antibiotic]["SS"],
                                                                    antibiotic2counts[antibiotic]["RS"],
                                                                    antibiotic2counts[antibiotic]["SR"],
                                                                    antibiotic2counts[antibiotic]["other"],
                                                                    s,
                                                                    sensitivity,
                                                                    specificity,
                                                                    precision))

    # sum all
    RR_sum = sum([antibiotic2counts[antibiotic]["RR"] for antibiotic in antibiotic2counts])
    SS_sum = sum([antibiotic2counts[antibiotic]["SS"] for antibiotic in antibiotic2counts])
    RS_sum = sum([antibiotic2counts[antibiotic]["RS"] for antibiotic in antibiotic2counts])
    SR_sum = sum([antibiotic2counts[antibiotic]["SR"] for antibiotic in antibiotic2counts])
    other_sum = sum([antibiotic2counts[antibiotic]["other"] for antibiotic in antibiotic2counts])

    sensitivity = round(RR_sum / (RR_sum + RS_sum), 3)
    specificity = round(SS_sum / (SS_sum + SR_sum), 3)
    precision = round(RR_sum / (RR_sum + SR_sum), 3)

    stats_f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("ALL",
                                                                RR_sum,
                                                                SS_sum,
                                                                RS_sum,
                                                                SR_sum,
                                                                other_sum,
                                                                s,
                                                                sensitivity,
                                                                specificity,
                                                                precision))






table_1 = snakemake.params["reference_table"] 
table_2 = snakemake.input["table_2"] 


disrepancies = snakemake.output["discrepancies"] 
FP = snakemake.output["FP"] 
FN = snakemake.output["FN"] 
counts = snakemake.output["counts"] 
stats = snakemake.output["stats"] 
   
COLUMNS = ["isoniazid","rifampicin","ethambutol","pyrazinamide"]

compare_tables(table_1, 
               table_2, 
               COLUMNS, 
               disrepancies,
               FP,
               FN,
               counts,
               stats)



