from Bio import SearchIO
import re 

ssearch_results = snakemake.input

o = open(snakemake.output[0], "w")
o.write("sample\tquery\thit\talignment_length\tpercent_identity\tevalue\tbitscore\n")
for ssearch in ssearch_results:
    print(ssearch)
    sample = re.search(".*/(.*)_ssearch.txt", ssearch).group(1)
    records = [i for i in SearchIO.parse(ssearch, 'fasta-m10')]
    for i,record in enumerate(records):
        for n,hit in enumerate(record.hits):
            
            description = hit.hsps[0].query_description
            query_start = hit.hsps[0].query_start
            query_end = hit.hsps[0].query_end
            hit_accession = hit.hsps[0].hit.id
            query_accession = hit.hsps[0].query.id
            alignment_length = hit.hsps[0].aln_span
            
            percent_identity = float(hit.hsps[0].ident_pct)
            evalue = hit.hsps[0].evalue
            bitscore = hit.hsps[0].bitscore
            # 16S_rRNA::NODE_52_length_5592_cov_1.520403:3723-5244(-)
            # 95266;tax=d:Bacteria,p:Bacteroidetes,c:Bacteroidia,o:Bacteroidales,f:Porphyromonadaceae,g:Parabacteroides,s:Parabacteroides
            
            o.write(f'{sample}\t{query_accession.split(":")[2]}\t{hit_accession.split("tax=")[1]}\t{alignment_length}\t{percent_identity}\t{evalue}\t{bitscore}\n')