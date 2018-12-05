
import pandas
import re
import numpy

# inputs
rgi_tsv_output = snakemake.input[0]
rgi_ontology = snakemake.input[1]
gene_depth_file = snakemake.input[2]
contig_gc_depth_file = snakemake.input[3]
samtools_depth = snakemake.input[4]
sample = snakemake.params[0]

# output
report_file = snakemake.output[0]

rgi_table = pandas.read_csv(rgi_tsv_output,
                            delimiter='\t',
                            header=0, index_col=0)

ontology_table = pandas.read_csv(rgi_ontology,
                                 delimiter='\t',
                                 header=0)

# parse rgi output file
contig2resistances = {}
with open(rgi_tsv_output, 'r') as f:
    for row in f:
        data = row.rstrip().split('\t')
        contig = '_'.join(data[1].split('_')[0:-1])
        if contig not in contig2resistances:
            contig2resistances[contig] = [re.sub("'", "", data[8])]
        else:
            contig2resistances[contig].append(re.sub("'", "", data[8]))

# get contig depth and GC
contig2gc_content = pandas.read_csv(contig_gc_depth_file,
                                    delimiter='\t',
                                    header=0).set_index("contig").to_dict()["gc_content"]

contig2median_depth = pandas.read_csv(contig_gc_depth_file,
                                      delimiter='\t',
                                      header=0).set_index("contig").to_dict()["median_depth"]

contig2size = pandas.read_csv(contig_gc_depth_file,
                              delimiter='\t',
                              header=0).set_index("contig").to_dict()["contig_size"]


# calculate rgi hit(s) sequencing depth based on position of the CDS in contigs
def parse_smatools_depth(samtools_depth):
    import pandas
    with open(samtools_depth, 'r') as f:
        table = pandas.read_csv(f, sep='\t', header=None, index_col=0)
    return table


def make_div(figure_or_data,
             include_plotlyjs=False,
             show_link=False,
             div_id=None):

    from plotly import offline

    div = offline.plot(figure_or_data,
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


def get_gc_coverage_data(contig2gc_content,
                         contig2resistances,
                         contig2size):
    import math

    bubble_data = {}
    bubble_data["with_resistance"] = {}
    bubble_data["with_resistance"]["hover_text"] = []
    bubble_data["with_resistance"]["contig_size"] = []
    bubble_data["with_resistance"]["GC"] = []
    bubble_data["with_resistance"]["depth"] = []

    bubble_data["without_resistance"] = {}
    bubble_data["without_resistance"]["hover_text"] = []
    bubble_data["without_resistance"]["contig_size"] = []
    bubble_data["without_resistance"]["GC"] = []
    bubble_data["without_resistance"]["depth"] = []

    for contig in contig2gc_content:
        depth = contig2median_depth[contig]
        gc = contig2gc_content[contig]
        if contig in contig2resistances:
            bubble_data["with_resistance"]["hover_text"].append("%s (%sbp): %s" % (contig,
                                                                                   contig2size[contig],
                                                                                   ','.join(contig2resistances[contig])))
            bubble_data["with_resistance"]["contig_size"].append(float(contig2size[contig]) + 10000)
            bubble_data["with_resistance"]["GC"].append(gc)
            bubble_data["with_resistance"]["depth"].append(depth)

        else:
            bubble_data["without_resistance"]["hover_text"].append("%s (%sbp)" % (contig,
                                                                                  contig2size[contig]))
            bubble_data["without_resistance"]["contig_size"].append(float(contig2size[contig]) + 10000)
            bubble_data["without_resistance"]["GC"].append(gc)
            bubble_data["without_resistance"]["depth"].append(depth)

    return bubble_data


def bubble_plot_gc_depth(bubble_data):
    import plotly.plotly as py
    import plotly.graph_objs as go

    import pandas as pd
    import math

    trace0 = go.Scatter(x=bubble_data["with_resistance"]["GC"],
                        y=bubble_data["with_resistance"]["depth"],
                        mode='markers',
                        name='With resistance(s)',
                        text=bubble_data["with_resistance"]["hover_text"],
                        marker=dict(symbol='circle',
                                    sizemode='area',
                                    sizeref=2.*max(bubble_data["without_resistance"]["contig_size"])/(40.**2),
                                    size=bubble_data["with_resistance"]["contig_size"],
                                    line=dict(width=2),
                                    )
                        )

    trace1 = go.Scatter(x=bubble_data["without_resistance"]["GC"],
                        y=bubble_data["without_resistance"]["depth"],
                        mode='markers',
                        name='Without resistance',
                        text=bubble_data["without_resistance"]["hover_text"],
                        marker=dict(symbol='circle',
                                    sizemode='area',
                                    sizeref=2.*max(bubble_data["without_resistance"]["contig_size"])/(40.**2),
                                    size=bubble_data["without_resistance"]["contig_size"],
                                    line=dict(width=2),
                                    )
                        )
    data = [trace0, trace1]
    layout = go.Layout(title='GC vs sequencing Depth plot',
                       width=1000,
                       height=700,
                       xaxis=dict(title='GC (%)',
                                  gridcolor='rgb(255, 255, 255)',
                                  # range=[2.003297660701705, 5.191505530708712],
                                  zerolinewidth=1,
                                  ticklen=5,
                                  gridwidth=2,
                                  ),
                       yaxis=dict(title='Depth',
                                  gridcolor='rgb(255, 255, 255)',
                                  # range=[36.12621671352166, 91.72921793264332],
                                  zerolinewidth=1,
                                  ticklen=5,
                                  gridwidth=2,
                                  ),
                       paper_bgcolor='rgb(243, 243, 243)',
                       plot_bgcolor='rgb(243, 243, 243)',
                       )
    fig = go.Figure(data=data, layout=layout)

    return (make_div(fig, div_id="bubble_plot"))


def resistance_table(rgi_table,
                     ontology_table,
                     samtools_depth_df):

    import numpy


    header = ["Contig",
              "ORF",
              "ARO",
              "Model Type",
              "Variant",
              "Cov (%)",
              "Identity (%)",
              "Score",
              "Score cutoff",
              "Family",
              "Mechanism",
              "Resistance",
              "Depth"]

    table_rows = []
    for n, one_resistance in rgi_table.iterrows():
        conting_str = one_resistance["Contig"].split('_')
        if len(conting_str) == 3:
            contig = '_'.join(conting_str[0:-1])
        elif len(conting_str) == 2:
            contig = one_resistance["Contig"]
        else:
            raise IOError("Unexpected contig format")

        orf_id = one_resistance["Contig"].split('_')[-1]
        gene_start = one_resistance["Start"]
        gene_end = one_resistance["Stop"]
        print(contig, n, one_resistance["Contig"])
        gene_depth = round(numpy.median(samtools_depth_df.loc[contig].iloc[gene_start:gene_end, 1]), 0)
        cutoff = one_resistance["Cut_Off"]
        name = one_resistance["Best_Hit_ARO"]
        identity = one_resistance["Best_Identities"]
        aro = one_resistance["ARO"]
        cov = one_resistance["Percentage Length of Reference Sequence"]
        bitscore = one_resistance["Best_Hit_Bitscore"]
        pass_bitscore = one_resistance["Pass_Bitscore"]
        mechanism = one_resistance["Resistance Mechanism"]
        snps = one_resistance["SNPs_in_Best_Hit_ARO"]
        model = one_resistance["Model_type"]
        family = one_resistance["AMR Gene Family"]
        antibio_res_list = list(ontology_table[ontology_table['Name'] == name]["Antibiotic resistance prediction"])
        antibio_res_class_list = list(ontology_table[ontology_table['Name'] == name]["Class"])

        resistance_code = ''
        for resistance, resistance_class in zip(antibio_res_list, antibio_res_class_list):
            # print(resistance, resistance_class)
            resistance_code += '%s (%s) </br>' % (resistance, resistance_class)

        table_rows.append([contig,
                           orf_id,
                           "%s (%s)" % (name, aro),
                           model,
                           snps,
                           cov,
                           identity,
                           bitscore,
                           pass_bitscore,
                           family,
                           mechanism,
                           resistance_code,
                           gene_depth])

    df = pandas.DataFrame(table_rows, columns=header)

    # cell content is truncated if colwidth not set to -1
    pandas.set_option('display.max_colwidth', -1)

    df_str = df.to_html(
        index=False,
        bold_rows=False,
        classes=["dataTable"],
        table_id="res_table",
        escape=False,
        border=0)

    return df_str.replace("\n", "\n" + 10 * " ")


STYLE = """
    <link rel="stylesheet" type="text/css" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css"/>
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.19/css/dataTables.bootstrap.min.css"/>
    <style type="text/css">
    body{font-family:Helvetica,arial,sans-serif;font-size:14px;line-height:1.6;padding-bottom:10px;background-color:#fff;color:#333;margin:0}body>div .section::before{content:"";display:block;height:80px;margin:-80px 0 0}#bubble-plot::before{margin:0}.topic-title{font-size:18pt}body>div>.section{margin-left:255px;margin-bottom:3em}div.section{margin-right:20px}#contents>p{display:none}button,li p.first{display:inline-block}#contents{margin-top:80px;padding-left:0;width:235px;background-color:#f1f1f1;height:100%;position:fixed;overflow:auto}#contents ul{list-style-type:none}#contents ul>li{font-size:14pt}#contents ul>li a:hover{color:#151d26}button,h1.title{color:#fff;background-color:#151d26}#contents ul>li>ul>li{font-size:12pt}h1.title{margin-top:0;position:fixed;z-index:10;padding:20px;width:100%}code,table tr:nth-child(2n),tt{background-color:#f8f8f8}.one-col{min-width:310px;height:500px;margin:0 auto}.two-col-left{height:300px;width:49%;float:left}.two-col-right{height:300px;width:49%;float:right}button{margin:0 5px 0 0;padding:5px 25px;font-size:18px;line-height:1.8;appearance:none;box-shadow:none;border-radius:3px;border:none}button:focus{outline:0}button:hover{background-color:#4183C4}button:active{background-color:#27496d}.legend-rect{width:20px;height:20px;margin-right:8px;margin-left:20px;float:left;-webkit-border-radius:2px;border-radius:2px}a{color:#4183C4;text-decoration:none}a.absent{color:#c00}a.anchor{padding-left:30px;margin-left:-30px;cursor:pointer;position:absolute;top:0;left:0;bottom:0}dl,dl dt,dl dt:first-child,hr,table,table tr{padding:0}table tr td,table tr th{border:1px solid #ccc;text-align:left;padding:6px 13px}h1,h2,h3,h4,h5,h6{margin:20px 0 10px;padding:0;font-weight:700;-webkit-font-smoothing:antialiased;cursor:text;position:relative}h1:hover a.anchor,h2:hover a.anchor,h3:hover a.anchor,h4:hover a.anchor,h5:hover a.anchor,h6:hover a.anchor{text-decoration:none}h1 code,h1 tt,h2 code,h2 tt,h3 code,h3 tt,h4 code,h4 tt,h5 code,h5 tt,h6 code,h6 tt{font-size:inherit}h1{font-size:28px;color:#151d26;border-bottom:1px solid #ccc}h2{font-size:24px;color:#000}h3{font-size:18px}h4{font-size:16px}dl dt,h5,h6{font-size:14px}h6{color:#777}blockquote,dl,li,ol,p,pre,table,ul{margin:15px 0}hr{background:url(http://tinyurl.com/bq5kskr) repeat-x;border:0;color:#ccc;height:4px}a:first-child h1,a:first-child h2,a:first-child h3,a:first-child h4,a:first-child h5,a:first-child h6{margin-top:0;padding-top:0}h1 p,h2 p,h3 p,h4 p,h5 p,h6 p{margin-top:0}dl dt{font-weight:700;font-style:italic;margin:15px 0 5px}blockquote>:first-child,dl dd>:first-child,dl dt>:first-child,table tr td :first-child,table tr th :first-child{margin-top:0}blockquote>:last-child,dl dd>:last-child,dl dt>:last-child{margin-bottom:0}dl dd{margin:0 0 15px;padding:0 15px}blockquote{border-left:4px solid #ddd;padding:0 15px;color:#777}table{border-spacing:0;border-collapse:collapse}table tr{border-top:1px solid #ccc;background-color:#fff;margin:0}table tr th{font-weight:700;margin:0}table tr td{margin:0}table tr td :last-child,table tr th :last-child{margin-bottom:0}img{max-width:100%}span.frame{display:block;overflow:hidden}span.frame>span{border:1px solid #ddd;display:block;float:left;overflow:hidden;margin:13px 0 0;padding:7px;width:auto}span.frame span img{display:block;float:left}span.frame span span{clear:both;color:#333;display:block;padding:5px 0 0}span.align-center{display:block;overflow:hidden;clear:both}span.align-center>span{display:block;overflow:hidden;margin:13px auto 0;text-align:center}span.align-center span img{margin:0 auto;text-align:center}span.align-right{display:block;overflow:hidden;clear:both}span.align-right>span{display:block;overflow:hidden;margin:13px 0 0;text-align:right}span.align-right span img{margin:0;text-align:right}span.float-left{display:block;margin-right:13px;overflow:hidden;float:left}span.float-left span{margin:13px 0 0}span.float-right{display:block;margin-left:13px;overflow:hidden;float:right}span.float-right>span{display:block;overflow:hidden;margin:13px auto 0;text-align:right}code,tt{margin:0 2px;padding:0 5px;white-space:nowrap;border:1px solid #eaeaea;border-radius:3px}pre code{margin:0;padding:0;white-space:pre;background:0 0}.highlight pre,pre{background-color:#f8f8f8;border:1px solid #ccc;font-size:13px;line-height:19px;overflow:auto;padding:6px 10px;border-radius:3px}pre code,pre tt{background-color:transparent;border:none}div#metadata{text-align:right}h1{line-height:1.6}.simple{padding-left:20px}.docutils.container{width:100%}
    .pull-left{
    .dataTables_filter {
       float: left !important;
    }
    mark.red {
        background-color: red;
        color: red;
    }
    mark.green {
        background-color: green;
        color: green;
    }
    </style>
    """

SCRIPT = """
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <script src="https://code.jquery.com/jquery-3.3.1.min.js"></script>
    <script src="https://cdn.datatables.net/1.10.19/js/jquery.dataTables.min.js"></script>
    <script src="https://cdn.datatables.net/1.10.19/js/dataTables.bootstrap.min.js"></script>

    <script>

    $(document).ready(function() {
        $('#res_table').DataTable( {
            dom: 'Bfrtip',
            "pageLength": 20,
            "searching": true,
            "bLengthChange": false,
            "paging":   true,
            "info": false,
            'rowCallback': function(row, data, index){
            if(data[6] < 90){
                             $(row).find('td:eq(6)').css('background-color', 'rgba(255, 0, 0, 0.6)');
                            }
            if(data[5] < 90){
                             $(row).find('td:eq(5)').css('background-color', 'rgba(0, 255, 0, 0.6)');
                            }
            },
        } );
    } );
    </script>

    """


def write_report(output_file,
                 STYLE,
                 SCRIPT,
                 table_rgi,
                 bubble_plot,
                 samtools_median_depth,
                 samtools_mean_depth):

    import io
    from docutils.core import publish_file, publish_parts
    from docutils.parsers.rst import directives

    if samtools_median_depth < 50:
        warning_type = 'danger'
    else:
        warning_type = 'info'

    report_str = f"""

.. raw:: html

    {SCRIPT}

    {STYLE}

=============================================================
RGI report
=============================================================

.. contents::
    :backlinks: none
    :depth: 2

Bubble plot
-----------


.. raw:: html

    {bubble_plot}


Table
------

.. raw:: html


    <div class="alert alert-{warning_type}" role="alert">
      Median depth: <strong> {samtools_median_depth} </strong> </br>
      Mean depth: <strong> {samtools_mean_depth} </strong>
    </div>

    <div class="alert alert-warning" role="alert">
      Genes with a hit <span class="label label-default">coverage &lt; 90%</span> are highlighted in <span class="label label-success">green</span> (if any) </br>
      Genes with an  <span class="label label-default">identity &lt; 90%</span> are highlighted in <span class="label label-danger">red</span> (if any)
    </div>

    {table_rgi}

"""
    with open(output_file, "w") as fh:
        publish_file(
            source=io.StringIO(report_str),
            destination=fh,
            writer_name="html",
            settings_overrides={"stylesheet_path": ""},
        )


plot_data = get_gc_coverage_data(contig2gc_content,
                                 contig2resistances,
                                 contig2size)

bubble_plot = bubble_plot_gc_depth(plot_data)

samtools_depth_df = parse_smatools_depth(samtools_depth)

samtools_mean_depth = round(numpy.mean(samtools_depth_df.iloc[:, 1]), 0)
samtools_median_depth = round(numpy.median(samtools_depth_df.iloc[:, 1]), 0)

rgi_table = resistance_table(rgi_table,
                             ontology_table,
                             samtools_depth_df)

write_report(report_file,
             STYLE,
             SCRIPT,
             rgi_table,
             bubble_plot,
             samtools_mean_depth,
             samtools_median_depth)
