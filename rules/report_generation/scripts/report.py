
import pandas
import re
from Bio import SeqIO
from plotly import offline


def checkm_table(checkm_table):
    
    df = pandas.read_csv(checkm_table, delimiter="\t", header=0)
    df = df[["Bin Id","Marker lineage","# markers", "Completeness","Contamination","Strain heterogeneity","0","1","2","3","4","5+"]]
    pandas.set_option('display.max_colwidth', -1)

    df_str = df.to_html(
        index=False,
        bold_rows=False,
        classes=["dataTable"],
        table_id="checkm_table",
        escape=False,
        border=0)

    return df_str.replace("\n", "\n" + 10 * " ")

def get_rrna_summary_table(raw_table, 
                           sample2scientific_name):
    
    df = pandas.read_csv(raw_table, delimiter="\t", header=0)
    
    # sample	query	hit	alignment_length	alignment_length	percent_identity	evalue	bitscore
    row_list = []
    sample2count = {}
    for n, row in df.iterrows():
        sample = str(row["sample"])
        contig = row["query"]
        BBH_taxnonomy = row["hit"].split(",")[1:]
        identity = row["percent_identity"]
        if sample not in sample2count:
            sample2count[sample] = 0
        sample2count[sample] += 1
        row_list.append([sample,
                         sample2count[sample], 
                         sample2scientific_name[sample], 
                         contig, 
                         BBH_taxnonomy, 
                         identity])
    
    header = ["sample", "N.", "expected", "contig", "BBH_taxnonomy", "identity"]
    df = pandas.DataFrame(row_list, columns=header)
    pandas.set_option('display.max_colwidth', -1)

    df_str = df.to_html(
        index=False,
        bold_rows=False,
        classes=["dataTable"],
        table_id="16S_table",
        escape=False,
        border=0)

    return df_str.replace("\n", "\n" + 10 * " ")


def get_centrifuge_table(centrifuge_links, sample2scientific_name):
    row_list = []
    for sample in centrifuge_links:
        sample_name = sample.split('/')[-2]
        table = pandas.read_csv(sample, delimiter="\t", names=["percent_root", "number_rooted", "number_assigned","rank", "taxid", "name"])

        total_alligned_reads = sum(table["number_assigned"])
        table["percent_assigned"] = round((table["number_assigned"] / total_alligned_reads) * 100, 2)

        #unclassified = table.iloc[0]
        hit_1 = table.iloc[0]
        hit_2 = table.iloc[1]
        hit_3 = table.iloc[2]

        row = [sample_name,
               sample2scientific_name[sample_name],
               "%s (%s)" % (hit_1["name"], hit_1["rank"]),
               hit_1["percent_assigned"],
               "%s (%s)" % (hit_2["name"], hit_2["rank"]),
               hit_2["percent_assigned"],
               "%s (%s)" % (hit_3["name"], hit_3["rank"]),
               hit_3["percent_assigned"],
               '<a href="%s">detail</a>' % '/'.join(sample.split('/')[1:])]

        row_list.append(row)
    header = ["sample", "expected", "taxon 1", "% reads", "taxon 2", "% reads", "taxon 3", "% reads", "detail"]
    df = pandas.DataFrame(row_list, columns=header)
    pandas.set_option('display.max_colwidth', -1)

    df_str = df.to_html(
        index=False,
        bold_rows=False,
        classes=["dataTable"],
        table_id="mash_table",
        escape=False,
        border=0)

    return df_str.replace("\n", "\n" + 10 * " ")

def get_mash_table(file_list, mash_detail, sample2scientific_name):
    '''
    sample
    hit 1: score (shared): name
    hit 2:
    hit 3:

    '''
    sample2detail_link = {}
    for link in mash_detail:
        sample = link.split("/")[-1].split('.')[0]
        sample2detail_link[str(sample)] = '/'.join(link.split("/")[1:])

    row_list = []
    for sample in file_list:
        sample_name = sample.split('/')[1]
        table = pandas.read_csv(sample, delimiter="\t", names=['score', 'shared_sketsh', 'e-value', 'description'])

        hit_1 = table.iloc[0]
        hit_2 = table.iloc[1]
        hit_3 = table.iloc[2]

        row = [sample_name,
               sample2scientific_name[sample_name],
               round(hit_1[0], 3),
               hit_1[1].split("/")[0],
               "%s..." % hit_1[3][0:25],
               round(hit_2[0], 3),
               hit_2[1].split("/")[0],
               "%s..." % hit_2[3][0:25],
               round(hit_3[0], 3),
               hit_3[1].split("/")[0],
               "%s..." % hit_3[3][0:25],
               '<a href="%s">detail</a>' % sample2detail_link[str(sample_name)]]

        row_list.append(row)
    header = ["sample", "expected", "score", "sketch", "description", "score", "sketch", "description","score", "sketch", "description", "detail"]
    df = pandas.DataFrame(row_list, columns=header)
    pandas.set_option('display.max_colwidth', -1)

    df_str = df.to_html(
        index=False,
        bold_rows=False,
        classes=["dataTable"],
        table_id="mash_table",
        escape=False,
        border=0)

    return df_str.replace("\n", "\n" + 10 * " ")

def get_multiqc_table(assembly_multiqc=False,
                      mapping_multiqc=False):

    mq_table = []
    if assembly_multiqc:
        mq_table.append(["MultiQC genome assemblie(s)", '<a href="%s">MultiQC</a>' % '/'.join(assembly_multiqc.split('/')[1:])])

    if mapping_multiqc:
        for n, multiqc in enumerate(mapping_multiqc):
            multiqc_link = '<a href="%s">MultiQC</a>' % '/'.join(multiqc.split('/')[1:])
            mq_table.append(["%s" % re.sub("_", " ", multiqc.split("/")[1]), multiqc_link])
    header = ["Name", "Link"]

    df = pandas.DataFrame(mq_table, columns=header)

    # cell content is truncated if colwidth not set to -1
    pandas.set_option('display.max_colwidth', -1)

    df_str = df.to_html(
        index=False,
        bold_rows=False,
        classes=["dataTable"],
        table_id="multiqc_table",
        escape=False,
        border=0)

    return df_str.replace("\n", "\n" + 10 * " ")


def get_core_genome_size(core_genome_bed):
    table = pandas.read_csv(core_genome_bed, delimiter="\t", names=['chromosome', 'start', 'end'])
    return sum(table["end"] - table["start"]+1)


def get_reference_genome_size(reference_genome_file):
    from Bio import SeqIO
    records_len = [len(i) for i in SeqIO.parse(open(reference_genome_file, 'r'), 'fasta')]
    return sum(records_len)


def quality_table(low_cov_fastas,
                  sample2gc,
                  sample2median_depth,
                  sampls2cumulated_size,
                  sampls2cumulated_size_filtered,
                  sample2n_contigs,
                  sample2scientific_name,
                  undetermined_snps_files=False,
                  core_genome_size=False,
                  low_cov_detail=False,
                  depth_cutoff=5):

    header = ["Strain id", "Scientific Name", "Contigs", f"Contigs depth <{depth_cutoff}", "GC", "Size >500pb (Mb)", f"Size >500bp & >{depth_cutoff}depth (Mb)", "Median Depth"]
    if low_cov_detail:
        sample2locov_path = {}
        for link in low_cov_detail:
            sample = link.split("/")[-1].split(".")[0]
            sample2locov_path[sample] = '/'.join(link.split("/")[1:])
    if undetermined_snps_files:
        # multiple files for each reference genome
        # first sort files
        reference2files = {}
        for one_file in undetermined_snps_files:
            # samples/{sample}/snps/{snp_caller}/{reference}/bwa/unknowns.tab
            reference_name = one_file.split('/')[4]
            if reference_name not in reference2files:
                reference2files[reference_name] = [one_file]
            else:
                reference2files[reference_name].append(one_file)
        reference2sample2n_unknown = {}
        for reference in reference2files:
            reference2sample2n_unknown[reference] = {}
            header.append("N. na SNP (%s)" % reference)
            if reference == 'cgMLST':
                header.append("Fraction core (%)")
            for f in reference2files[reference]:
                sample = f.split('/')[1]
                try:
                    reference2sample2n_unknown[reference][sample] = len(pandas.read_csv(f, delimiter='\t'))
                except pandas.errors.EmptyDataError:
                    reference2sample2n_unknown[reference][sample] = 0
    cov_table = []
    for fasta in low_cov_fastas:

        sample = fasta.split("/")[1]
        try:
            with open(fasta, 'r') as f:
                n_records = len([i for i in SeqIO.parse(f, 'fasta')])
        except ValueError:
            n_records = 0
        if n_records > 0:
            low_cov_str = '%s (<a href="%s">detail<a>)' % (n_records, sample2locov_path[sample])
        else:
            low_cov_str = n_records

        if undetermined_snps_files:

            tmp_lst = [sample,
                       sample2scientific_name[sample],
                       sample2n_contigs[sample],
                       low_cov_str,
                       sample2gc[sample],
                       round(sampls2cumulated_size[sample] / 1000000, 2),
                       round(sampls2cumulated_size_filtered[sample] / 1000000, 2),
                       sample2median_depth[sample]]
            for reference in reference2files:
                tmp_lst.append(reference2sample2n_unknown[reference][sample])
                if reference == "cgMLST":
                    tmp_lst.append(round((float(reference2sample2n_unknown[reference][sample]) / core_genome_size) * 100, 2))

        else:
            tmp_lst = [sample,
                       sample2scientific_name[sample],
                       sample2n_contigs[sample],
                       low_cov_str,
                       sample2gc[sample],
                       round(sampls2cumulated_size[sample] / 1000000, 2),
                       round(sampls2cumulated_size_filtered[sample] / 1000000, 2),
                       sample2median_depth[sample]]
        cov_table.append(tmp_lst)

    if len(cov_table) > 0:
        df = pandas.DataFrame(cov_table, columns=header)

        # cell content is truncated if colwidth not set to -1
        pandas.set_option('display.max_colwidth', -1)

        df_str = df.to_html(
            index=False,
            bold_rows=False,
            classes=["dataTable"],
            table_id="cov_table",
            escape=False,
            border=0)

        return df_str.replace("\n", "\n" + 10 * " ")
    else:
        return 'No sample with low coverage contigs'


def qualimap_table(qualimap_links, self_mapping=False):

    cov_table = []
    sample2ref2link = {}
    for qualimap in qualimap_links:
        search = re.search('report/qualimap/(.*)/.*/(.*)/qualimapReport.html', qualimap)
        sample = search.group(1)
        ref = search.group(2)
        if sample in ref:
            ref = "self mapping"
        mod_path = '/'.join(qualimap.split('/')[1:])
        if sample not in sample2ref2link:
            sample2ref2link[sample] = {}
            sample2ref2link[sample][ref] = "<a href=%s>Qualimap report</a>" % mod_path
        else:
            sample2ref2link[sample][ref] = "<a href=%s>Qualimap report</a>" % mod_path
    for sample in sample2ref2link:
        cov_table.append([sample] + [sample2ref2link[sample][i] for i in sample2ref2link[sample]])

    header = ["Strain id"] + ["Ref: %s" % i for i in sample2ref2link[sample]]
    df = pandas.DataFrame(cov_table, columns=header)

    # cell content is truncated if colwidth not set to -1
    pandas.set_option('display.max_colwidth', -1)

    df_str = df.to_html(
        index=False,
        bold_rows=False,
        classes=["dataTable"],
        table_id="qualimap_table",
        escape=False,
        border=0)

    return df_str.replace("\n", "\n" + 10 * " ")


def virulence_table(virulence_reports,
                    blast_files,
                    ordered_samples,
                    fasta_files=False):

    header = ["Strain id","Number of VFs","Virulence Report"]

    sample2n_VFs = {}
    for n, sample in enumerate(ordered_samples):
        if not fasta_files:
            sample2n_VFs[sample] = len(blast_files[n])
        else:
            sample2n_VFs[sample] = len([i for i in SeqIO.parse(open(blast_files[n], 'r'), "fasta")])

    vf_data = []
    # report/virulence/VFDB/{sample}_report.html
    report_template = '<a href="virulence/%s/%s_report.html">Details</a>'
    for report in virulence_reports:
        sample = re.search('report/virulence/.*/(.*)_report.html', report).group(1)
        type = re.search('report/virulence/(.*)/.*_report.html', report).group(1)
        vf_data.append([sample,
                        sample2n_VFs[sample],
                        report_template % (type, sample)])

    df = pandas.DataFrame(vf_data, columns=header)

    # cell content is truncated if colwidth not set to -1
    pandas.set_option('display.max_colwidth', -1)

    df_str = df.to_html(
        index=False,
        bold_rows=False,
        classes=["dataTable"],
        table_id="VF_table",
        escape=False,
        border=0)

    return df_str.replace("\n", "\n" + 10 * " ")


def resistance_table(resistance_reports):

    header = ["Strain id", "Resistance Report"]

    rgi_data = []
    report_template = '<a href="resistance/%s_rgi_report.html">RGI report</a>'
    for report in resistance_reports:
        sample = re.search('report/resistance/(.*)_rgi_report.html', report).group(1)
        rgi_data.append([sample,
                      report_template % sample])

    df = pandas.DataFrame(rgi_data, columns=header)

    # cell content is truncated if colwidth not set to -1
    pandas.set_option('display.max_colwidth', -1)

    df_str = df.to_html(
        index=False,
        bold_rows=False,
        classes=["dataTable"],
        table_id="RGI_table",
        escape=False,
        border=0)

    return df_str.replace("\n", "\n" + 10 * " ")


def plot_minimum_spanning_tree():
    pass

def make_div(figure_or_data, include_plotlyjs=False, show_link=False, div_id=None):
    div = offline.plot(
        figure_or_data,
        include_plotlyjs=include_plotlyjs,
        show_link=show_link,
        output_type="div",
    )
    if ".then(function ()" in div:
        div = f"""{div.partition(".then(function ()")[0]}</script>"""
    if div_id:
        import re

        try:
            existing_id = re.findall(r'id="(.*?)"|$', div)[0]
            div = div.replace(existing_id, div_id)
        except IndexError:
            pass
    return div


def plot_heatmap_snps(mat, id):
    import pandas
    import scipy.cluster.hierarchy as hc
    import plotly.figure_factory as ff
    import plotly.graph_objs as go
    import numpy

    m = pandas.read_csv(mat, delimiter='\t', header=0, index_col=0)

    link = hc.linkage(m.values, method='centroid')
    o1 = hc.leaves_list(link)
    mat = m.iloc[o1, :]
    mat2 = mat.iloc[:, o1[::-1]].iloc[::-1]

    nodes = mat2.index.tolist()

    data = ff.create_annotated_heatmap(
        z= mat2.values,  # squareform(m.values)
        x=['S_' + str(i) for i in mat2.columns.values.tolist()],
        y=['S_' + str(i) for i in mat2.index.tolist()],  # need to reverse order for y axis
        colorscale='Reds'
    )

    layout = go.Layout(
        title='SNP heatmap'
    )

    fig = go.Figure(data=data, layout=layout)
    fig.layout.margin.update({"l": 20 + (max([len(i) for i in nodes]) * 7),
                              "r": 0,
                              "b": 20,
                              "t": 60,
                              "pad": 10,
                              })
    return make_div(fig, div_id=id)


def link_list2dico(link_list, label, index_sample, add_cgMLST=False):
    reference2sample2link = {}
    for link in link_list:
        data = link.split("/")
        reference = data[2]
        sample = data[index_sample].split(".")[0]
        if reference not in reference2sample2link:
            reference2sample2link[reference] = {}
        reference2sample2link[reference][sample] = '<a href="%s">%s</a>' % ('/'.join(data[1:]), label)
    if add_cgMLST and not 'cgMLST' in reference2sample2link:
        reference2sample2link["cgMLST"] = {}
        for sample in reference2sample2link[reference]:
            reference2sample2link["cgMLST"][sample] = "-"
    return reference2sample2link


def get_snp_detail_table(snp_link_list, indel_link_list):
    '''
    report/snps/TATRas-control_assembled_genome/bwa/gatk_gvcfs/TATRas-control.html
    report/snps/TATRas-control_assembled_genome/bwa/gatk_gvcfs/TATRas-mutant-A.html
    report/snps/TATRas-control_assembled_genome/bwa/gatk_gvcfs/TATRas-mutant-B.html
    report/snps/cgMLST/bwa/gatk_gvcfs/TATRas-control.html
    report/snps/cgMLST/bwa/gatk_gvcfs/TATRas-mutant-A.html
    report/snps/cgMLST/bwa/gatk_gvcfs/TATRas-mutant-B.html
    '''

    reference2sample2snp_link = link_list2dico(snp_link_list, "snps", 5)
    if len(indel_link_list) > 0:
        reference2sample2indel_link = link_list2dico(indel_link_list, "del", 4, add_cgMLST=True)

    reference_list = list(reference2sample2snp_link.keys())
    header = ["Sample"] + ["Reference: %s" % i for i in reference_list]
    rows = []
    for sample in list(reference2sample2snp_link[reference_list[0]].keys()):
        try:
            rows.append([sample] + ["%s / %s" % (reference2sample2snp_link[ref][sample], reference2sample2indel_link[ref][sample]) for ref in reference_list])
        except:
            rows.append([sample] + ["%s / %s" % (reference2sample2snp_link[ref][sample], "-") for ref in reference_list])
    df = pandas.DataFrame(rows, columns=header)

    # cell content is truncated if colwidth not set to -1
    pandas.set_option('display.max_colwidth', -1)

    df_str = df.to_html(
        index=False,
        bold_rows=False,
        classes=["dataTable"],
        table_id="snp_detail_table",
        escape=False,
        border=0)

    return df_str.replace("\n", "\n" + 10 * " ")
