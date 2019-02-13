
import pandas

report = snakemake.input["report"]

report_table = pandas.read_csv(report, names=["percent_root", "number_root", "number_assigned", "rank", "taxid", "name"], delimiter='\t')

output = open(snakemake.output["formatted_out"], 'w')

total_alligned_reads = sum(report_table["number_assigned"])
report_table["percent_assigned"] = round((report_table["number_assigned"] / total_alligned_reads) * 100, 2)
report_table = report_table.sort_values(by=['percent_assigned'], ascending=False)

report_table.to_csv(output, sep="\t")
