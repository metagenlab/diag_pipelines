import pandas
import openpyxl
import re


out_xlsx = snakemake.output[0]

excels = [pandas.ExcelFile(name) for name in snakemake.input]

writer = pandas.ExcelWriter(out_xlsx)


for j in excels:
    for sheet in j.sheet_names:
        df = j.parse(sheet)
        df.to_excel(writer, sheet, index=True)
        writer.save()
