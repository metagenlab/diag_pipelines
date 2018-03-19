from Bio import SeqIO

def extract_translation(filename, locus, out, missing):
    with open(out, "w") as outfile:
        for record in SeqIO.parse(filename, "genbank"):
            version=(record.id)
            for f in record.features:
                if f.type=="gene":
                    if "old_locus_tag" in f.qualifiers:
                        if locus & set(f.qualifiers["old_locus_tag"]):
                            outfile.write(version+"\t"+str(f.location.start)+"\t"+str(f.location.end)+"\n")
                            locus.remove(list(locus & set(f.qualifiers["old_locus_tag"]))[0])
    open(missing, "w").write("\n".join(list(locus)))

gbk = snakemake.input["gbk"]
locus_tags = [x.strip() for x in list(open(snakemake.input["locus_list"]))]
extract_translation(gbk, set(locus_tags), snakemake.output["bed"], snakemake.output["problematic"])
