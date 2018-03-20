from Bio import SeqIO

def extract_translation(filename, locus, out):
    with open(out, "w") as outfile:
        for record in SeqIO.parse(filename, "genbank"):
            version=(record.id)
            for f in record.features:
                if f.type=="gene":
                    if locus & set(f.qualifiers["locus_tag"]):
                        loc = locus & set(f.qualifiers["locus_tag"])
                        if f.location.strand == -1:
                            outfile.write(">"+str(list(loc)[0])+"\n"+str(record.seq[f.location.start:f.location.end].reverse_complement())+"\n")
                        else:
                            outfile.write(">"+str(list(loc)[0])+"\n"+str(record.seq[f.location.start:f.location.end])+"\n")

    
gbk = snakemake.input["gbk"]
locus_tags = [x.strip() for x in list(open(snakemake.input["locus_list"]))]
print(locus_tags)
extract_translation(gbk, set(locus_tags), snakemake.output["genes"])
