
import pandas
import re
from Bio import SeqIO
from plotly import offline


def get_core_genome_size(core_genome_bed):
    table = pandas.read_csv(core_genome_bed, delimiter="\t", names=['chromosome', 'start', 'end'])
    return sum(table["end"]-table["start"]+1)

def get_reference_genome_size(reference_genome_file):
    from Bio import SeqIO
    records_len = [len(i) for i in SeqIO.parse(open(reference_genome_file, 'r'), 'fasta')]
    return sum(records_len)


def coverage_table(low_cov_fastas):

    header = ["Strain id","Number of contigs"]

    cov_table = []
    for fasta in low_cov_fastas:
        sample = re.search('samples/(.*)/assembly/spades/coverage_filtered/contigs_500bp_low_coverage.fasta', fasta).group(1)
        try:
            with open(fasta, 'r') as f:
                n_records = len(SeqIO.read(f, 'fasta'))
        except ValueError:
            continue
        cov_table.append([sample, n_records])

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

def qualimap_table(qualimap_links):

    header = ["Strain id","Number of contigs"]

    cov_table = []
    for qualimap in qualimap_links:
        sample = re.search('report/qualimap/(.*)/bwa/.*_assembled_genome/qualimapReport.html', qualimap).group(1)
        mod_path = '/'.join(qualimap.split('/')[1:])
        cov_table.append([sample, "<a href=%s>Qualimap report</a>" % mod_path])

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
                    ordered_samples):

    header = ["Strain id","Number of VFs","Virulence Report"]

    sample2n_VFs = {}
    for n, sample in enumerate(ordered_samples):
        sample2n_VFs[sample] = len(blast_files[n])

    vf_data = []
    report_template = '<a href="virulence/%s_VFDB_report.html">VFDB report</a>'
    for report in virulence_reports:
        sample = re.search('report/virulence/(.*)_VFDB_report.html', report).group(1)
        vf_data.append([sample,
                      sample2n_VFs[sample],
                      report_template % sample])

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

    header = ["Strain id","Resistance Report"]

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


def plot_heatmap_snps(mat):
    import pandas
    import scipy.cluster.hierarchy as hc
    import plotly.figure_factory as ff
    import plotly.graph_objs as go

    m = pandas.read_csv(mat, delimiter='\t', header=0, index_col=0)
    link = hc.linkage(m.values, method='centroid')
    o1 = hc.leaves_list(link)

    mat = m.iloc[o1, :]
    mat = mat.iloc[:, o1[::-1]]

    nodes = ['S_' + str(i) for i in mat.index]

    data = ff.create_annotated_heatmap(
        z=mat.values,  # squareform(m.values)
        x=nodes,
        y=nodes,
        colorscale='Reds'
    )

    layout = go.Layout(
        title='SNP heatmap'
    )

    fig = go.Figure(data=data, layout=layout)

    return make_div(fig, div_id="heatmapPlot")
