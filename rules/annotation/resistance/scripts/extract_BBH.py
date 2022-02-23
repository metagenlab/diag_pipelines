
from Bio import SeqIO
import pandas 

gbk = SeqIO.parse(snakemake.input[1], "genbank")

ssearch_result = pandas.read_csv(snakemake.input[0], header = None, names = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"], sep="\t")

# contig	gene	start	end	depth	ratio_assembly	contig_depth	contig_ratio_depth	contig_length
gene_depth_file = pandas.read_csv(snakemake.input["gene_depth_file"], sep="\t", header=0, index_col="gene")

fasta_db = SeqIO.parse(snakemake.params["db_fasta"], "fasta")
id_cutoff = float(snakemake.params["id_cutoff"])

# >gi|40950501|gb|AAR97884.1|AFA-1| class A beta-lactamase
# >WP_063645611.1|AFA-2| class A beta-lactamase
accession2data = {record.id: [record.description, record.seq] for record in fasta_db}

locus2data = {}
for record in gbk:
    for feature in record.features:
        if feature.type == 'CDS':

            locus2data[feature.qualifiers["locus_tag"][0]] = {"start": int(feature.location.start), 
                                                              "stop": int(feature.location.end), 
                                                              "strand": feature.location.strand, 
                                                              "nucl_seq": str(feature.extract(record.seq)),
                                                              "aa_seq": str(feature.qualifiers["translation"][0]),
                                                              "contig": record.name, 
                                                              "depth": int(gene_depth_file.loc[feature.qualifiers["locus_tag"][0],"depth"]), 
                                                              "contig_depth": int(gene_depth_file.loc[feature.qualifiers["locus_tag"][0],"contig_depth"]),
                                                              "contig_length": int(gene_depth_file.loc[feature.qualifiers["locus_tag"][0],"contig_length"])}





header = {
"qseqid": "ORF_ID",		
"bitscore": "Bitscore",
"pident": "Percent_identity",
"ID": "ID",
"Model_ID": "Model_ID",	
"Nudged": "Nudged",
"Note": "Note"	
}

ssearch_result = ssearch_result.rename(columns=header)

ssearch_result["reference_db"] = "BLDB"
ssearch_result["Sample"] = snakemake.params[0]
ssearch_result["Contig"] = None
ssearch_result["Cut_Off"] = None
ssearch_result["Pass_Bitscore"] = None
ssearch_result["Best_Hit"] = None
ssearch_result["ARO"] = None
ssearch_result["Model_type"] = "protein homolog model"
ssearch_result["SNPs"] = None
ssearch_result["Drug_class"] = "Beta-Lactam"
ssearch_result["Mechanism"] = None
ssearch_result["Note"] = None
ssearch_result["Cut_Off"] = None

locus_list = set()
keep = []
for n, row in ssearch_result.iterrows():
    # keep only BBH
    if row["ORF_ID"] not in locus_list:
        # filter low identity hits
        if float(row.Percent_identity) < id_cutoff:
            continue
        locus_list.add(row["ORF_ID"])

        hit_label = row["sseqid"]
        hit_name = hit_label.split("|")[-2]
        hit_accession = hit_label.split("|")[-3]
        hit_description = accession2data[row["sseqid"]][0].split("|")[-1].strip()
        hit_seq_length = len(accession2data[row["sseqid"]][1])

        orf_id = row.ORF_ID.split(" ")[0]

        row["ORF_ID"] = orf_id
        row["Best_hit"] = hit_name
        row["AMR_family"] = hit_description
        align_length = int(row.qend)-int(row.qstart)
        row["Percent_coverage"] = round((align_length/float(hit_seq_length))*100,2)
        row["Start"] = locus2data[orf_id]["start"]
        row["Contig"] = locus2data[orf_id]["contig"]
        row["Stop"] = locus2data[orf_id]["stop"]
        row["Strand"] = locus2data[orf_id]["strand"]
        row["Protein_sequence"] = locus2data[orf_id]["aa_seq"]
        row["Nucleotide_sequence"] = locus2data[orf_id]["nucl_seq"]
        row["Gene_depth"] = locus2data[orf_id]["depth"]
        row["Contig_depth"] = locus2data[orf_id]["contig_depth"]
        row["Contig_length"] = locus2data[orf_id]["contig_length"]

        if row["Percent_coverage"] > 90 and row["Percent_identity"] == 100:
            row["Cut_Off"] = 'Perfect'
        else:
            row["Cut_Off"] = 'Strict'

        keep.append(row)

ssearch_result_filtered = pandas.DataFrame(keep)
print("SSEARCH OUT", ssearch_result_filtered)
if len(ssearch_result_filtered) > 0:
    ssearch_result_filtered["Start"] = ssearch_result_filtered["Start"].astype(int)
    ssearch_result_filtered["Stop"] = ssearch_result_filtered["Stop"].astype(int)
    ssearch_result_filtered["Strand"] = ssearch_result_filtered["Strand"].astype(int)
    ssearch_result_filtered["Gene_depth"] = ssearch_result_filtered["Gene_depth"].astype(int)
    ssearch_result_filtered["Contig_depth"] = ssearch_result_filtered["Contig_depth"].astype(int)
    ssearch_result_filtered["Contig_length"] = ssearch_result_filtered["Contig_length"].astype(int)
    ssearch_result_filtered[["Sample", "reference_db", "ORF_ID", "Contig", "Contig_length", "Contig_depth", "Start", "Stop", "Strand", "Gene_depth","Best_hit", "ARO", "Bitscore","Pass_Bitscore", "Cut_Off" ,"Percent_identity", "Percent_coverage", "AMR_family", "Model_type", "Mechanism", "SNPs", "Drug_class", "Protein_sequence", "Nucleotide_sequence", "Note"]].to_csv(snakemake.output[0], sep="\t", index=None)
else:
    # empty dataframe
    with open(snakemake.output[0], "w") as f:
        f.write('\t'.join(["Sample", "reference_db", "ORF_ID", "Contig", "Contig_length", "Contig_depth", "Start", "Stop", "Strand", "Gene_depth","Best_hit", "ARO", "Bitscore","Pass_Bitscore", "Cut_Off" ,"Percent_identity", "Percent_coverage", "AMR_family", "Model_type", "Mechanism", "SNPs", "Drug_class", "Protein_sequence", "Nucleotide_sequence", "Note"]))
