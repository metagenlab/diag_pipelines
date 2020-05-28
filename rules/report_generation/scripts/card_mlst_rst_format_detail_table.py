#!/usr/bin/env python
    
def parse_rgi(rgi_file_list,
              output_file,
              query_cov_cutoff=50):
    import pandas
    import re

    o = open(output_file, 'w')

    header = [
        "Sample",
        "Contig_id",
        "Start",
        "Stop",
        "Strand",
        "Best_hit",
        "Cut_Off",
        "Percent_identity",
        "Percent_coverage",
        "AMR_family",
        "Model_type",
        "Mechanism",
        "SNPs",
        "Drug_class",
        "Note",
        "Protein_sequence",
        "Nucleotide_sequence"
    ]
    o.write("\t".join(header)+'\n')

    for rgi_file in rgi_file_list:
        sample = rgi_file.split("/")[1]
        t = pandas.read_csv(rgi_file, sep="\t", header=0)
        # ORF_ID	
        # Contig	
        # Start	
        # Stop	
        # Orientation	
        # Cut_Off	
        # Pass_Bitscore	
        # Best_Hit_Bitscore	
        # Best_Hit_ARO	
        # Best_Identities	
        # ARO	
        # Model_type	
        # SNPs_in_Best_Hit_ARO	
        # Other_SNPs	
        # Drug Class	
        # Resistance Mechanism	
        # AMR Gene Family	
        # Predicted_DNA	
        # Predicted_Protein	
        # CARD_Protein_Sequence	
        # Percentage Length of Reference Sequence	
        # ID	
        # Model_ID	
        # Nudged	
        # Note


        for n, row in t.iterrows():
            contig_prefix, contig_number, orf_number = row["Contig"].split("_")
            contig_id = f'{contig_prefix}_{contig_number}'
            gene = row["Best_Hit_ARO"]
            model_type = row["Model_type"]
            mechanism = row["Resistance Mechanism"]
            identity = row["Best_Identities"]
            snps = row["SNPs_in_Best_Hit_ARO"]
            drug_class = row["Drug Class"]
            family = row["AMR Gene Family"]
            query_cov = row["Percentage Length of Reference Sequence"]
            note = row["Note"]
            protein_seq = row["Predicted_Protein"]
            dna_seq = row["Predicted_DNA"]
            start = row["Start"]
            end = row["Stop"]
            strand = row["Orientation"]
            cutoff = row["Cut_Off"]
            nudged = row["Nudged"]
            note = row["Note"]
            
            if float(query_cov) < query_cov_cutoff:
                print("Skipping low cov entry: %s (%s %%)" % (gene, query_cov))
                continue

            o.write(f'{sample}\t{contig_id}\t{start}\t{end}\t{strand}\t{gene}\t{cutoff}' \
                  f'\t{identity}\t{query_cov}\t{family}\t{model_type}\t{mechanism}' \
                  f'\t{snps}\t{drug_class}\t{note}\t{protein_seq}\t{dna_seq}\n')
            

rgi_tables = snakemake.input["rgi_files"]
output_file = snakemake.output[0]

sample2rgi = parse_rgi(rgi_tables,output_file)
