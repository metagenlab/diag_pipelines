import pandas
import collections
import itertools


def remove_redundancy_from_blastp(results):
    protein_duplicates = [item for item, count in collections.Counter(list(results["matching_sequence"])).items() if count > 1]
    if len(protein_duplicates)==0:
        return(results)
    else:
        for i in protein_duplicates:
            matching_rows = results.loc[data["matching_sequence"]==i]
            results = results.drop(matching_rows.index)
            results = results.append(matching_rows.loc[matching_rows["percentage_identity"].idxmax()])
        return(results)


def remove_redundancy_from_tblastn(results):
    contig_duplicates = [item for item, count in collections.Counter(list(results["matching_sequence"])).items() if count > 1]
    if len(contig_duplicates)==0:
        return(results)
    else:
        for i in contig_duplicates:
            to_be_kept=[]
            to_be_deleted=[]
            overlapping_hits = {}
            matching_rows = results.loc[data["matching_sequence"]==i]
            for i in itertools.combinations(matching_rows.index, 2):
                ranges = [range(sorted((int(x[0]),int(x[1])))[0], sorted((int(x[0]),int(x[1])))[1]) for x in matching_rows.loc[i, ["alignment_start_on_matching_sequence", "aligment_end_on_matching_sequence"]].values]
                if set(ranges[0]) & set(ranges[1]):
                    to_be_kept.append(matching_rows.loc[i, "percentage_identity"].idxmax())
                    to_be_deleted.append(matching_rows.loc[i, "percentage_identity"].idxmin())
            if not set(to_be_kept) & set(to_be_deleted):
                results = results.drop(to_be_deleted)
            else:
                raise ValueError("Problem during removing overlapping results from the tblastn results, one hit overlapped with at least two other hits!")
        return(results)



results ={}
data = pandas.read_csv(snakemake.input[0], sep="\t", header=0, index_col=0)
print(data)

blastp_matches = data.loc[data["algorithm"]=="blastp"]

tblastn_matches = data.loc[data["algorithm"]=="tblastn"]

blastp_matches = remove_redundancy_from_blastp(blastp_matches)
tblastn_matches = remove_redundancy_from_tblastn(tblastn_matches)

tblastn_matches.append(blastp_matches).to_csv(snakemake.output[0], sep="\t")
