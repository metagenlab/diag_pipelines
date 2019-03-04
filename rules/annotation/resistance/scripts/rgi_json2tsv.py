

import json

rgi_tsv_output = snakemake.output[0]
rgi_json = snakemake.input[0]

"""
ORF_ID	contig_1_83 # 83178 # 83840 # 1 # ID=1_83;partial=00;start_type=TTG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.626 ==> contig
Contig	contig_1_83 ==> 'orf_from'
Start	83178 ==> 'orf_start'
Stop	83840 ==> 'orf_end'
Orientation	+ ==> 'orf_strand'
Cut_Off	Strict ==> 'type_match'
Pass_Bitscore 420 ==> 'pass_bitscore'
Best_Hit_Bitscore 435.261 ==> 'bit_score'
Best_Hit_ARO	mtrA ==> 'ARO_name' / 'model_name'
Best_Identities	99.55 ==> 'perc_identity'
ARO	3000816 ==> 'ARO_accession'
Model_type	protein homolog model ==>'model_type'
SNPs_in_Best_Hit_ARO n/a ==>
Other_SNPs n/a ==> !!!!!!!!!!!!!!!!!
Drug Class macrolide antibiotic; penam ==> RM
Resistance Mechanism antibiotic efflux ==> RM
AMR Gene Family resistance-nodulation-cell division (RND) antibiotic efflux pump ==> RM
Predicted_DNA ==> 'orf_dna_sequence'
Predicted_Protein ==> 'orf_prot_sequence'
CARD_Protein_Sequence ==> sequence_from_db'
Percentage Length of Reference Sequence 96.49 ==> round(('match'/'sequence_from_broadstreet')*100, 2)
ID gnl|BL_ORD_ID|1463|hsp_num:0  ==> hit
Model_ID 1559 ==> 'model_id',

snp	{'original': 'L', 'change': 'R', 'position': 511}

snp = ['snp']["original"] + ['snp']["position"] + ['snp']["change"]

[, , , , ,, 'model_type_id', ,  'pass_evalue', , , , 'ARO_category', 'evalue', , 'max_identities', 'cvterm_id', 'query','match' , ', 'sequence_from_broadstreet', 'dna_sequence_from_broadstreet', 'query_start', 'query_end', , , ]
"""

o = open(rgi_tsv_output, 'w')

header = ["ORF_ID",
          "Contig",
          "Start",
          "Stop",
          "Orientation",
          "Cut_Off",
          "Pass_Bitscore",
          "Best_Hit_Bitscore",
          "Best_Hit_ARO",
          "Best_Identities",
          "ARO",
          "Model_type",
          "SNPs_in_Best_Hit_ARO",
          "Other_SNPs",
          "Drug Class",
          "Resistance Mechanism",
          "AMR Gene Family",
          "Predicted_DNA",
          "Predicted_Protein",
          "CARD_Protein_Sequence",
          "Percentage Length of Reference Sequence",
          "ID",
          "Model_ID"]

o.write("\t".join(header)+'\n')
with open(rgi_json) as f:
    data = json.loads(f.read())
    gene2snps = {}
    for contig in data.keys():
        for hit in data[contig]:
            if not 'hsp_num:0' in hit:
                continue
            amr_data = {}
            if data[contig][hit]["type_match"] != "Loose":
                for key in data[contig][hit]["ARO_category"]:
                    if data[contig][hit]["ARO_category"][key]["category_aro_class_name"] == 'Resistance Mechanism':
                        amr_data["Resistance Mechanism"] = data[contig][hit]["ARO_category"][key]['category_aro_name']
                    elif data[contig][hit]["ARO_category"][key]["category_aro_class_name"] == 'AMR Gene Family':
                        amr_data["AMR Gene Family"] = data[contig][hit]["ARO_category"][key]["category_aro_name"]
                    else:
                        #print("-------- ARO: %s ---------" % data[contig][hit]['ARO_accession'])
                        #print ("%s\t%s\t%s\t%s" % (contig, hit, key, data[contig][hit]["ARO_category"][key]))
                        amr_data["Resistance Mechanism"] = data[contig][hit]["ARO_category"][key]['category_aro_name']

                amr_data["ORF_ID"] = contig
                amr_data["Contig"] = data[contig][hit]['orf_from']
                amr_data["Start"] = data[contig][hit]['orf_start']
                amr_data["Stop"] = data[contig][hit]['orf_end']
                amr_data["Orientation"] = data[contig][hit]['orf_strand']
                amr_data["Cut_Off"] = data[contig][hit]['type_match']
                amr_data["Pass_Bitscore"] = data[contig][hit]['pass_bitscore']
                amr_data["Best_Hit_Bitscore"] = data[contig][hit]['bit_score']
                amr_data["Best_Hit_ARO"] = data[contig][hit]['ARO_name']
                amr_data["Best_Identities"] = data[contig][hit]['perc_identity']
                amr_data["ARO"] = data[contig][hit]['ARO_accession']
                amr_data["Model_type"] = data[contig][hit]['model_type']
                amr_data["SNPs_in_Best_Hit_ARO"] = "n/a"
                amr_data["Other_SNPs"] = "n/a"
                amr_data["Drug Class"] = "n/a"
                amr_data["Predicted_DNA"] = data[contig][hit]['orf_dna_sequence']
                amr_data["Predicted_Protein"] = data[contig][hit]['orf_prot_sequence']
                amr_data["CARD_Protein_Sequence"] = data[contig][hit]['sequence_from_db']
                seq1 = data[contig][hit]['match']
                seq2 = data[contig][hit]['sequence_from_broadstreet']
                #print(amr_data["Best_Hit_ARO"])
                #print ("seq1", seq1)
                #print("seq2", seq2)
                #print(data[contig][hit])
                try:
                    amr_data["Percentage Length of Reference Sequence"] = round((len(seq1)/len(seq2))*100, 2)
                except:
                    amr_data["Percentage Length of Reference Sequence"] = 0
                amr_data["ID"] = hit
                amr_data["Model_ID"] = data[contig][hit]['model_id']

                if not "snp" in data[contig][hit]:
                    row = [amr_data[i] for i in header]
                    o.write('\t'.join([str(i) for i in row])+'\n')
                else:
                    new_snp = data[contig][hit]['snp']["original"] + str(data[contig][hit]['snp']["position"]) + data[contig][hit]['snp']["change"]
                    if not amr_data["Contig"] in gene2snps:
                        amr_data["SNPs_in_Best_Hit_ARO"] = [new_snp]
                        gene2snps[contig] = amr_data
                        gene2snps[contig]["Predicted_Protein"] = "n/a"
                        gene2snps[contig]["Predicted_DNA"] = "n/a"
                        gene2snps[contig]["CARD_Protein_Sequence"] = "n/a"
                    else:
                        gene2snps[contig]["SNPs_in_Best_Hit_ARO"].append(new_snp)
    for gene in gene2snps:
        amr_data = gene2snps[gene]
        amr_data["SNPs_in_Best_Hit_ARO"]= ','.join(amr_data["SNPs_in_Best_Hit_ARO"])
        row = [amr_data[i] for i in header]
        o.write('\t'.join([str(i) for i in row])+'\n')
