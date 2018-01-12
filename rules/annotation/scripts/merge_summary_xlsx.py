import pandas
import openpyxl
import re

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]


snakemake.input.sort(key=natural_keys)

out_xlsx = snakemake.output[0]

excels = [pandas.ExcelFile(name) for name in snakemake.input]

writer = pandas.ExcelWriter(out_xlsx)

for j in excels:
    for sheet in j.sheet_names:
        df = j.parse(sheet)
        df.to_excel(writer, sheet, index=False)
        writer.save()
