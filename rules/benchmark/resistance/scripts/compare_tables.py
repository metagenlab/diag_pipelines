#!/usr/bin/env python

import pandas

def compare_tables(t1,t2, columns, rgi):

    RR_count = 0
    SS_count = 0
    RS_count = 0 
    SR_count = 0
    other = 0

    table_a = pandas.read_csv(t1, index_col=0, sep="\t")
    table_b = pandas.read_csv(t2, index_col=0, sep="\t")
    print("file1\tfile2\tsample\tdrug\tsuceptibility1\ttsuceptibility2")
    drug2counts = {} 
    for column in columns:
        drug2counts[column] = {"reference" : { "R" : 0, "S": 0, "other": 0}  , "tool" :  { "R" : 0, "S": 0, "other": 0}}  

    for sample in table_a.index:
        for column in columns:

            sa = table_a.at[sample, column]
            colref = column
            if rgi == True:
                if column == "rifampicin":
                    column = "rifamycin"
                elif column == "ethambutol":
                    column = "polyamine"

            
            sb = table_b.at[sample, column]

            if sa == 'R':
                drug2counts[colref]["reference"]["R"] +=1 
            elif sa == 'S':
                drug2counts[colref]["reference"]["S"] +=1
            elif sa not in ["S", "R"] :
                drug2counts[colref]["reference"]["other"] +=1
            if sb == 'R':
                drug2counts[colref]["tool"]["R"] +=1 
            elif sb == 'S':
                drug2counts[colref]["tool"]["S"] +=1
            elif sb not in ["S", "R"] :
                drug2counts[colref]["tool"]["other"] +=1


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
    for drug in drug2counts:
        
        print("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (drug,
                                              drug2counts[drug]["reference"]["R"],
                                              drug2counts[drug]["reference"]["S"],
                                              drug2counts[drug]["reference"]["other"],
                                              drug2counts[drug]["tool"]["R"],
                                              drug2counts[drug]["tool"]["S"],
                                              drug2counts[drug]["tool"]["other"]))

    print("\n OVERVIEW:")
    print(f"RR\t{RR_count}")
    print(f"SS\t{SS_count}")
    print(f"RS\t{RS_count}")
    print(f"SR\t{SR_count}")
    print(f"Other\t{other}")


table_1 = snakemake.params["reference_table"] 
table_2 = snakemake.input["table_2"] 


# RGI: rifampicin = rifamycin
# RGI: ethambutol = polyamine  
if "rgi" in table_2:
    rgi = True
else:
    rgi = False
    
COLUMNS = ["isoniazid","rifampicin","ethambutol","pyrazinamide"]

compare_tables(table_1, table_2, COLUMNS, rgi=rgi)



