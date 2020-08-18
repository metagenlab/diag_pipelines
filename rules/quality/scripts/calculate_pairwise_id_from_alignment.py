
from Bio import AlignIO
import numpy as np 

def pairewise_identity(seq1, seq2):
    "return identity calculated as: n identical sites/n aligned sites (gaps non included)"
    A=list(seq1.upper())
    B=list(seq2.upper())
    identical_sites = 0
    aligned_sites = 0
    gaps=0
    for n in range(0, len(A)):
        if str(A[n]) != "-" and str(B[n]) != "-":
            aligned_sites+=1
        else:
            continue
        if str(A[n])==str(B[n]):
            identical_sites+=1
    try:
        identity = 100*(identical_sites/float(aligned_sites))
    except ZeroDivisionError:
        identity = 0

    #print 'identical_sites', identical_sites
    #print 'aligned_sites', aligned_sites

    if aligned_sites == 0:
        # return false if align length = 0
        return False, False
    else:
        return 100*(identical_sites/float(aligned_sites)), aligned_sites


def get_identity_matrix_from_multiple_alignment(alignment):
    # [len(alignment)+1, len(alignment)]
    identity_matrix = np.chararray((len(alignment)+1, len(alignment)+1), itemsize=30)
    matrix_aligned_sites = np.chararray((len(alignment)+1, len(alignment)+1), itemsize=30)
    identity_matrix[0,0] = "-"
    # first column = locus tags
    for x in range(0,len(alignment)):
        identity_matrix[x+1, 0] = alignment[x].name
        identity_matrix[0, x+1] = alignment[x].name

        for y in range(x, len(alignment)):
            identity, n_aligned_sites = pairewise_identity(alignment[x], alignment[y])
            #print "identity", identity
            identity_matrix[y+1, x+1] = round(identity, 2)
            identity_matrix[x+1, y+1] = round(identity, 2)
            matrix_aligned_sites[x+1, y+1] = n_aligned_sites
            matrix_aligned_sites[x+1, y+1] = n_aligned_sites
    #print identity_matrix
    return identity_matrix, matrix_aligned_sites


alignment = snakemake.input[0]
o = open(snakemake.output[0], "w")


align = AlignIO.read(alignment, "fasta")

m, aligned_sites = get_identity_matrix_from_multiple_alignment(align)

for i in range(1, len(m[0,:])-1):
    for y in range(i+1, len(m[0,:])):
        if i==y:
            continue
        if m[0,i].decode("utf-8").split('|')[0] == m[0,y].decode("utf-8").split('|')[0]:
            o.write("%s\t%s\t%s\t%s" %('_'.join(m[0,i].decode("utf-8").split('|')[0:2]),
                                    '_'.join(m[0,y].decode("utf-8").split('|')[0:2]),
                                    m[i,y].decode("utf-8"),
                                    aligned_sites[i,y].decode("utf-8") ))