
import venn
import pandas

rgi = snakemake.input["rgi"]
mykrobe = snakemake.input["mykrobe"]
tbprofiler = snakemake.input["tbprofiler"]
walker = snakemake.input["walker"]

antibio = snakemake.params["antibio"]

o = open(snakemake.output[0], 'w')

# rgi	ERR2512416	macrolide	Mycobacterium avium 23S rRNA with mutation conferring resistance to clarithromycin	2274	A2274G	SNP	-
# custom	ERR2512419	rifampicin	rpoB	450	S450L	SNP	R
header = [
    "tool",
    "sample",
    "antibiotic",
    "gene",
    "position",
    "variant",
    "variant_type",
    "phenotype"
]


rgi_FP_counts = pandas.read_csv(rgi, sep="\t", header=None, names=header).query('antibiotic == "%s" & phenotype == "S"' % antibio)[["gene", "variant"]].groupby(["gene","variant"]).size().to_dict()
mykrobe_FP_counts = pandas.read_csv(mykrobe, sep="\t", header=None, names=header).query('antibiotic == "%s" & phenotype == "S"' % antibio)[["gene", "variant"]].groupby(["gene","variant"]).size().to_dict()
tbprofiler_FP_counts = pandas.read_csv(tbprofiler, sep="\t", header=None, names=header).query('antibiotic == "%s" & phenotype == "S"' % antibio)[["gene", "variant"]].groupby(["gene","variant"]).size().to_dict()
walker_FP_counts = pandas.read_csv(walker, sep="\t", header=None, names=header).query('antibiotic == "%s" & phenotype == "S"' % antibio)[["gene", "variant"]].groupby(["gene","variant"]).size().to_dict()


nr_key_list = list(set(list(rgi_FP_counts.keys()) + list(mykrobe_FP_counts.keys()) + list(tbprofiler_FP_counts.keys())+ list(walker_FP_counts.keys())))

o.write("gene\tvariant\trgi_freq\tmykrobe_freq\ttbprofiler_freq\twalker_freq\n")
for one_var in nr_key_list:
    if one_var in rgi_FP_counts:
        rgi_freq = rgi_FP_counts[one_var]
    else:
        rgi_freq = 0
    
    if one_var in mykrobe_FP_counts:
        mykrobe_freq = mykrobe_FP_counts[one_var]
    else:
        mykrobe_freq = 0

    if one_var in tbprofiler_FP_counts:
        tbprofiler_freq = tbprofiler_FP_counts[one_var]
    else:
        tbprofiler_freq = 0
        
    if one_var in walker_FP_counts:
        walker_freq = walker_FP_counts[one_var]
    else:
        walker_freq = 0
    o.write(f"{one_var[0]}\t{one_var[1]}\t{rgi_freq}\t{mykrobe_freq}\t{tbprofiler_freq}\t{walker_freq}\n")