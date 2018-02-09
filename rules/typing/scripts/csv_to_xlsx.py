import pandas
from io import StringIO

writer = pandas.ExcelWriter(snakemake.output[0])

df = pandas.read_csv(snakemake.input[0], sep="\t", header=None, index_col=0)
header = pandas.Series(df.drop([1,2], axis=1).iloc[0]).str.replace("\(", "").str.replace("\)", "").str.replace("\d", "")
values = df.drop([1,2], axis=1).replace(to_replace="[a-z]", value="", regex=True).replace(to_replace="[A-Z]", value="", regex=True).replace(to_replace="\(", value="", regex=True).replace(to_replace="\)", value="", regex=True)
final = pandas.concat([df[1], df[2], values], axis=1)
final.to_excel(writer, header=["MLST schema", "ST types"] + list(header), sheet_name="mlst summary")
writer.save()
