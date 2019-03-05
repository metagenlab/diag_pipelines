
import pandas

report = snakemake.input["report"]
detail = snakemake.input["detail"]


detail_table = pandas.read_csv(detail, header=0, delimiter='\t')
report_table = pandas.read_csv(report, header=0, delimiter='\t')

output = open(snakemake.output["formatted_out"], 'w')

total_reads = len(detail_table["readID"])
print("total",len(detail_table["readID"]))
print("nr",len(set(detail_table["readID"])))

unclassified = len(detail_table["seqID"][detail_table["seqID"]=="unclassified"])

unclassified_proportion = round((unclassified/total_reads)*100, 2)

# name	taxID	taxRank	genomeSize	numReads numUniqueReads	abundance
report_table["proportion"] = round((report_table["numReads"] / total_reads) * 100, 2)
report_table = report_table.sort_values(by=['proportion'], ascending=False)

output.write("name\ttaxid\trank\tn_reads\tpercent_reads\n")
output.write("unclassified\t-\t-\t%s\t%s\n" % (unclassified, unclassified_proportion))
for n, row in report_table.iterrows():
    name = row["name"]
    taxid = row["taxID"]
    rank = row["taxRank"]
    n_reads = row["numReads"]
    percent_reads = round((row["numReads"] / total_reads) * 100, 2)
    output.write("%s\t%s\t%s\t%s\t%s\n" % (name, taxid, rank, n_reads, percent_reads))
