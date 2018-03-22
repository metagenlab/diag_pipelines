import pandas
import re
import Bio

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

df = pandas.read_csv(snakemake.input[0], sep="\t")


rgi_res = pandas.read_csv(snakemake.input[1], sep="\t")

#If we only used RGI, we report all the protein associated with resistance from the RGI results
if all([x=="rgi" for x in list(df["Software"])]):
    with open(snakemake.output[1], "w") as f_p, open(snakemake.output[0], "w") as f_dna:
        written = []
        for i in df.index:
            id_seq = snakemake.wildcards["sample"] + "_" + df.loc[i, "Gene"]
            if id_seq not in written:
                match = rgi_res.loc[rgi_res["Best_Hit_ARO"]==df.loc[i, "Gene"].replace("_", " "), "Predicted_DNA"]
                dna_seq = Seq(str(match.values[0]))
                dna_record = SeqRecord(dna_seq, id_seq, description="")
                SeqIO.write(dna_record, f_dna, "fasta")
                prot_seq = rgi_res.loc[rgi_res["Best_Hit_ARO"]==df.loc[i, "Gene"].replace("_", " "), "Predicted_Protein"].values[0]
                #Careful here, I don't know how RGI results look if there is no predicted protein (rrs)
                if prot_seq:
                    prot_seq = Seq("M" + prot_seq[1:])
                    prot_record = SeqRecord(prot_seq, id=id_seq, description="")
                    SeqIO.write(prot_record, f_p, "fasta")
            written.append(id_seq)



#If we used RGI and mykrobe, we report the protein that were found by both softwares
else:
    with open(snakemake.output[1], "w") as f_p, open(snakemake.output[0], "w") as f_dna:
        for i in df.loc[df['Software'] == "mykrobe", "Gene"]:
            #the matching between the results of rgi and mykrobe is done here, as each software names gene differently, this can be highly unreliable
            #a complete mapping between mykrobe and rgi gene names should be performed...
            match=set([string for string in df.loc[df['Software'] == "rgi", "Gene"] if re.search(i, string.replace("(",""))])
            written = []
            for j in match:
                id_seq = snakemake.wildcards["sample"] + "_" + j
                if id_seq not in written:
                    dna_seq = Seq(rgi_res.loc[rgi_res["Best_Hit_ARO"]==j.replace("_", " "), "Predicted_DNA"].values[0], generic_dna)
                    dna_record= SeqRecord(dna_seq, id=id_seq, description="")
                    SeqIO.write(dna_record, f_dna, "fasta")
                    prot_seq = rgi_res.loc[rgi_res["Best_Hit_ARO"]==j.replace("_", " "), "Predicted_Protein"].values[0]
                    #Careful here, I don't know how RGI results look if there is no predicted protein (rrs)
                    if prot_seq:
                        prot_seq = Seq("M" + prot_seq[1:])
                        prot_record = SeqRecord(prot_seq, id=id_seq, description="")
                        SeqIO.write(prot_record, f_p, "fasta")
                    written.append(id_seq)
                        


