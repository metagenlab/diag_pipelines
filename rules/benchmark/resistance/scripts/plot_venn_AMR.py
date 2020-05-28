
import venn
import pandas

rgi = snakemake.input["rgi"]
mykrobe = snakemake.input["mykrobe"]
tbprofiler = snakemake.input["tbprofiler"]
walker = snakemake.input["walker"]

antibio = snakemake.params["antibio"]

if antibio == "all":
    rgi_set = set(pandas.read_csv(rgi, sep="\t", header=None, names=["db","sample","antibiotic"])["sample"].tolist())
    mykrobe_set = set(pandas.read_csv(mykrobe, sep="\t", header=None, names=["db","sample","antibiotic"])["sample"].tolist())
    tbprofiler_set = set(pandas.read_csv(tbprofiler, sep="\t", header=None, names=["db","sample","antibiotic"])["sample"].tolist())
    walker_set = set(pandas.read_csv(walker, sep="\t", header=None, names=["db","sample","antibiotic"])["sample"].tolist())
else:
    rgi_set = set(pandas.read_csv(rgi, sep="\t", header=None, names=["db","sample","antibiotic"]).query('antibiotic == "%s"' % antibio)["sample"].tolist())
    mykrobe_set = set(pandas.read_csv(mykrobe, sep="\t", header=None, names=["db","sample","antibiotic"]).query('antibiotic == "%s"' % antibio)["sample"].tolist())
    tbprofiler_set = set(pandas.read_csv(tbprofiler, sep="\t", header=None, names=["db","sample","antibiotic"]).query('antibiotic == "%s"' % antibio)["sample"].tolist())
    walker_set = set(pandas.read_csv(walker, sep="\t", header=None, names=["db","sample","antibiotic"]).query('antibiotic == "%s"' % antibio)["sample"].tolist())


labels = venn.get_labels([rgi_set, mykrobe_set, tbprofiler_set, walker_set], fill=['number'])

fig, ax = venn.venn4(labels, names=['rgi', 'mykrobe', 'tbprofiler', 'walker'])

fig.savefig(snakemake.output[0])