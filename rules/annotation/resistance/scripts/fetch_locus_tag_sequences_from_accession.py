from Bio import SeqIO

def extract_translation(filename, locus, genes_out, shift):
    with open(genes_out, "w") as genes:
        for record in SeqIO.parse(filename, "genbank"):
            version=(record.id)
            for f in record.features:
                if f.type=="gene":
                    if locus & set(f.qualifiers["locus_tag"]):
                        loc = locus & set(f.qualifiers["locus_tag"])
                        if f.location.strand == -1:
                            genes.write(">"+str(list(loc)[0])+"\n"+str(record.seq[f.location.start-shift:f.location.end+shift].reverse_complement())+"\n")
                        else:
                            genes.write(">"+str(list(loc)[0])+"\n"+str(record.seq[f.location.start-shift:f.location.end+shift])+"\n")

    
gbk = snakemake.input["gbk"]
locus_tags = [x.strip() for x in list(open(snakemake.input["locus_list"]))]
extract_translation(gbk, set(locus_tags), snakemake.output["genes"], snakemake.params["upstream_downstream_size"])

