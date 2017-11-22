import Bio
import Bio.SeqIO
import Bio.AlignIO

def read_alignment(filename):
    return(Bio.AlignIO.parse(filename, "fasta"))

def return_seqRec(aliRec, cor):
    seqs={}
    for i in aliRec:
        for j in i:
            seqs[cor[j.id]]=j.seq
    return(seqs)


def return_number_of_snps(seqRec):
    s = "\t"+"\t".join(seqRec.keys())+"\n"
    for i in seqRec.keys():
        s += i + "\t"
        for j in seqRec.keys():
            nogaps = sum(aa1 != "-" and aa2 != "-" for aa1, aa2 in zip(seqRec[i], seqRec[j]))
            matches = sum(aa1 == aa2 and aa1 != "-" and aa1 != "n" and aa1 != "N" for aa1, aa2 in zip(seqRec[i], seqRec[j]))
            s += "{:d} \t".format(nogaps-matches)
        s += "\n"
    return(s)

corres = snakemake.params.corres
corres["Reference"] = snakemake.wildcards.ref

with open(snakemake.output[0], "w") as myfile:
    myfile.write(return_number_of_snps(return_seqRec(read_alignment(snakemake.input[0]), corres)))
