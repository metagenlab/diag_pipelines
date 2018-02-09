import pandas
import collections
import itertools

def remove_redundancy_from_blastp(results):
    #get duplicate entries from the list of matching proteins
    protein_duplicates = [item for item, count in collections.Counter(list(results["matching_sequence"])).items() if count > 1]
    if len(protein_duplicates)==0:
        return(results)
    else:
        for i in protein_duplicates:
            matching_rows = results.loc[data["matching_sequence"]==i]
            #remove all the duplicated entries
            results = results.drop(matching_rows.index)
            #put back only the best entry based on the percentage identity
            results = results.append(matching_rows.loc[matching_rows["percentage_identity"].idxmax()])
        return(results)


def remove_redundancy_from_tblastn(results):
    #get duplicate entries from the list of matching contigs
    contig_duplicates = [item for item, count in collections.Counter(list(results["matching_sequence"])).items() if count > 1]
    if len(contig_duplicates)==0:
        return(results)
    else:
        for i in contig_duplicates:
            to_be_kept = []
            to_be_deleted = []
            matching_rows = results.loc[data["matching_sequence"]==i]
            # analyse every pair of hits that are on the same contig 
            for i in itertools.combinations(matching_rows.index, 2):
                # get list with the ranges of the two hits
                list_start_end = [ sorted((int(x[0]),int(x[1]))) for x in matching_rows.loc[i, ["alignment_start_on_matching_sequence", "aligment_end_on_matching_sequence"]].values ]
                ranges = [ ranges(x[0], x[1]) for x in list_start_end ]
                # if the ranges overlap, we mark the best hit for keeping and the worst for deleting
                if set(ranges[0]) & set(ranges[1]):
                    to_be_kept.append(matching_rows.loc[i, "percentage_identity"].idxmax())
                    to_be_deleted.append(matching_rows.loc[i, "percentage_identity"].idxmin())
            # if a hit intersects with two other hits we can have problem, if it was marked for deletion during one comparison and marked for keeping in the other
            # this case will result in failure, should be taken care of if it ever happens
            if not set(to_be_kept) & set(to_be_deleted):
                results = results.drop(to_be_deleted)
            else:
                raise ValueError("Problem during removing overlapping results from the tblastn results, one hit overlapped with at least two other hits!")
        return(results)


data = pandas.read_csv(snakemake.input[0], sep="\t", header=0, index_col=0)

blastp_matches = data.loc[data["algorithm"]=="blastp"]

tblastn_matches = data.loc[data["algorithm"]=="tblastn"]

blastp_matches = remove_redundancy_from_blastp(blastp_matches)
tblastn_matches = remove_redundancy_from_tblastn(tblastn_matches)

tblastn_matches.append(blastp_matches).to_csv(snakemake.output[0], sep="\t")
