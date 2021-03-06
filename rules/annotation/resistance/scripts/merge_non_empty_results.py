import pandas
import openpyxl
import re

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split('(\d+)', text) ]


snakemake.input["xlsx"].sort(key=natural_keys)

out_xlsx = snakemake.output[0]

excels = [pandas.ExcelFile(name) for name in snakemake.input["xlsx"]]

writer = pandas.ExcelWriter(out_xlsx)

for j in excels:
    for sheet in j.sheet_names:
        df = j.parse(sheet)
        if df.shape[0]:
            df.to_excel(writer, sheet_name=sheet, index=True)
try:            
    writer.save()
except IndexError:
    pandas.DataFrame(["Empty"]).to_excel(writer)
    writer.save()

