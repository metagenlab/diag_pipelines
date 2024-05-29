import pandas
from report import quality_table, resistance_table
import report
import io
from docutils.core import publish_file
from Bio import SeqIO

multiqc_assembly = snakemake.input["multiqc_assembly"]
rgi_overview = '/'.join(snakemake.input["rgi_overview"].split('/')[1:])
ordered_samples = snakemake.params["samples"]
resistance_reports = snakemake.input["resistance_reports"]
low_cov_fastas = snakemake.input["low_cov_fastas"]
output_file = snakemake.output[0]
high_cov_fastas = snakemake.input["high_cov_fastas"]
high_cov_fastgs = snakemake.input["high_cov_fastgs"]
contig_gc_depth_file_list = snakemake.input["contig_gc_depth_file_list"]
mash_results = snakemake.input["mash_results"]  # ok
qualimap_reports = snakemake.input["qualimap_reports"]
low_cov_detail = snakemake.input["low_cov_detail"]
mash_detail = snakemake.input["mash_detail"]
sample2scientific_name = snakemake.params["sample_table"].to_dict()["ScientificName"]
centrifuge_tables = snakemake.input["centrifuge_tables"]
checkm_table = snakemake.input["checkm_table"]
rrna_classification_file = snakemake.input["rrna_classification_file"]
depth_cutoff = snakemake.params["depth_cutoff"]

STYLE = """
    <link rel="stylesheet" type="text/css" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css"/>
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.19/css/dataTables.bootstrap.min.css"/>
    <style type="text/css">
    body{font-family:Helvetica,arial,sans-serif;font-size:14px;line-height:1.6;padding-bottom:10px;background-color:#fff;color:#333;margin:0}body>div .section::before{content:"";display:block;height:80px;margin:-80px 0 0}#quality-control::before{margin:0}.topic-title{font-size:18pt}body>div>.section{margin-left:255px;margin-bottom:3em}div.section{margin-right:20px}#contents>p{display:none}button,li p.first{display:inline-block}#contents{margin-top:80px;padding-left:0;width:235px;background-color:#f1f1f1;height:100%;position:fixed;overflow:auto}#contents ul{list-style-type:none}#contents ul>li{font-size:14pt}#contents ul>li a:hover{color:#151d26}button,h1.title{color:#fff;background-color:#151d26}#contents ul>li>ul>li{font-size:12pt}h1.title{margin-top:0;position:fixed;z-index:10;padding:20px;width:100%}code,table tr:nth-child(2n),tt{background-color:#f8f8f8}.one-col{min-width:310px;height:500px;margin:0 auto}.two-col-left{height:300px;width:49%;float:left}.two-col-right{height:300px;width:49%;float:right}button{margin:0 5px 0 0;padding:5px 25px;font-size:18px;line-height:1.8;appearance:none;box-shadow:none;border-radius:3px;border:none}button:focus{outline:0}button:hover{background-color:#4183C4}button:active{background-color:#27496d}.legend-rect{width:20px;height:20px;margin-right:8px;margin-left:20px;float:left;-webkit-border-radius:2px;border-radius:2px}a{color:#4183C4;text-decoration:none}a.absent{color:#c00}a.anchor{padding-left:30px;margin-left:-30px;cursor:pointer;position:absolute;top:0;left:0;bottom:0}dl,dl dt,dl dt:first-child,hr,table,table tr{padding:0}table tr td,table tr th{border:1px solid #ccc;text-align:left;padding:6px 13px}h1,h2,h3,h4,h5,h6{margin:20px 0 10px;padding:0;font-weight:700;-webkit-font-smoothing:antialiased;cursor:text;position:relative}h1:hover a.anchor,h2:hover a.anchor,h3:hover a.anchor,h4:hover a.anchor,h5:hover a.anchor,h6:hover a.anchor{text-decoration:none}h1 code,h1 tt,h2 code,h2 tt,h3 code,h3 tt,h4 code,h4 tt,h5 code,h5 tt,h6 code,h6 tt{font-size:inherit}h1{font-size:28px;color:#151d26;border-bottom:1px solid #ccc}h2{font-size:24px;color:#000}h3{font-size:18px}h4{font-size:16px}dl dt,h5,h6{font-size:14px}h6{color:#777}blockquote,dl,li,ol,p,pre,table,ul{margin:15px 0}hr{background:url(http://tinyurl.com/bq5kskr) repeat-x;border:0;color:#ccc;height:4px}a:first-child h1,a:first-child h2,a:first-child h3,a:first-child h4,a:first-child h5,a:first-child h6{margin-top:0;padding-top:0}h1 p,h2 p,h3 p,h4 p,h5 p,h6 p{margin-top:0}dl dt{font-weight:700;font-style:italic;margin:15px 0 5px}blockquote>:first-child,dl dd>:first-child,dl dt>:first-child,table tr td :first-child,table tr th :first-child{margin-top:0}blockquote>:last-child,dl dd>:last-child,dl dt>:last-child{margin-bottom:0}dl dd{margin:0 0 15px;padding:0 15px}blockquote{border-left:4px solid #ddd;padding:0 15px;color:#777}table{border-spacing:0;border-collapse:collapse}table tr{border-top:1px solid #ccc;background-color:#fff;margin:0}table tr th{font-weight:700;margin:0}table tr td{margin:0}table tr td :last-child,table tr th :last-child{margin-bottom:0}img{max-width:100%}span.frame{display:block;overflow:hidden}span.frame>span{border:1px solid #ddd;display:block;float:left;overflow:hidden;margin:13px 0 0;padding:7px;width:auto}span.frame span img{display:block;float:left}span.frame span span{clear:both;color:#333;display:block;padding:5px 0 0}span.align-center{display:block;overflow:hidden;clear:both}span.align-center>span{display:block;overflow:hidden;margin:13px auto 0;text-align:center}span.align-center span img{margin:0 auto;text-align:center}span.align-right{display:block;overflow:hidden;clear:both}span.align-right>span{display:block;overflow:hidden;margin:13px 0 0;text-align:right}span.align-right span img{margin:0;text-align:right}span.float-left{display:block;margin-right:13px;overflow:hidden;float:left}span.float-left span{margin:13px 0 0}span.float-right{display:block;margin-left:13px;overflow:hidden;float:right}span.float-right>span{display:block;overflow:hidden;margin:13px auto 0;text-align:right}code,tt{margin:0 2px;padding:0 5px;white-space:nowrap;border:1px solid #eaeaea;border-radius:3px}pre code{margin:0;padding:0;white-space:pre;background:0 0}.highlight pre,pre{background-color:#f8f8f8;border:1px solid #ccc;font-size:13px;line-height:19px;overflow:auto;padding:6px 10px;border-radius:3px}pre code,pre tt{background-color:transparent;border:none}div#metadata{text-align:right}h1{line-height:1.6}.simple{padding-left:20px}.docutils.container{width:100%}
    .pull-left{
    .dataTables_filter {
       float: left !important;
    }
    </style>
    """

SCRIPT = """
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <script src="https://code.jquery.com/jquery-3.3.1.min.js"></script>
    <script src="https://cdn.datatables.net/1.10.19/js/jquery.dataTables.min.js"></script>
    <script src="https://cdn.datatables.net/1.10.19/js/dataTables.bootstrap.min.js"></script>
    <script src="https://unpkg.com/cytoscape/dist/cytoscape.min.js"></script>
    <script src="https://unpkg.com/webcola/WebCola/cola.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/cytoscape-cola@2.2.4/cytoscape-cola.min.js"></script>

    <script>

    $(document).ready(function() {
        $('#RGI_table').DataTable( {
            dom: 'Bfrtip',
            "pageLength": 20,
            "searching": true,
            "bLengthChange": false,
            "paging":   true,
            "info": false
        } );
    } );

    $(document).ready(function() {
        $('#cov_table').DataTable( {
            dom: 'Bfrtip',
            "pageLength": 20,
            "searching": true,
            "bLengthChange": false,
            "paging":   true,
            "info": false
        } );
    } );
    $(document).ready(function() {
        $('#mash_table').DataTable( {
            dom: 'Bfrtip',
            "pageLength": 20,
            "searching": true,
            "bLengthChange": false,
            "paging":   true,
            "info": false,
            'rowCallback': function(row, data, index){
                $(row).find('td:eq(2)').css('background-color', 'rgba(255, 0, 0, 0.2)');
                $(row).find('td:eq(3)').css('background-color', 'rgba(255, 0, 0, 0.2)');
                $(row).find('td:eq(4)').css('background-color', 'rgba(255, 0, 0, 0.2)');
                $(row).find('td:eq(5)').css('background-color', 'rgba(0,128,0, 0.2)');
                $(row).find('td:eq(6)').css('background-color', 'rgba(0,128,0, 0.2)');
                $(row).find('td:eq(7)').css('background-color', 'rgba(0,128,0, 0.2)');
                $(row).find('td:eq(8)').css('background-color', 'rgba(128,128,128, 0.2)');
                $(row).find('td:eq(9)').css('background-color', 'rgba(128,128,128, 0.2)');
                $(row).find('td:eq(10)').css('background-color', 'rgba(128,128,128, 0.2)');
            },
        } );
    } );
    </script>

    """

# get contig depth and GC
sample2gc = {}
sample2median_depth = {}
sampls2cumulated_size = {}
sample2n_contigs = {}
sampls2cumulated_size_filtered = {}
sample2no_n_contigs = {}

for one_table in contig_gc_depth_file_list:
    table = pandas.read_csv(one_table,
                            delimiter='\t',
                            header=0,
                            index_col=0)

    data_whole_gnome = table.loc["TOTAL"]
    n_contigs = len(table["gc_content"])-1

    # samples/5965/quality/mapping/bwa/5965_assembled_genome/contig_gc_depth_500bp_high_coverage.tab
    sample = one_table.split('/')[1]

    sample2gc[sample] = data_whole_gnome["gc_content"]
    sample2median_depth[sample] = data_whole_gnome["mean_depth"]
    sampls2cumulated_size[sample] = data_whole_gnome["contig_size"]
    sampls2cumulated_size_filtered[sample] = int(table.query('median_depth>=5 & contig != "TOTAL"')[["contig_size"]].sum())
    sample2n_contigs[sample] = n_contigs

#Counting contigs with no neightbor in graph (higher perc of dead ends in assembly graph means poorer assembly)
for fastgs in high_cov_fastgs:
    sample = fastgs.split("/")[1]
    tot = 0
    no_neighbor_contigs = 0
    for header in SeqIO.parse(fastgs, "fasta"):
        if not ":EDGE" in header.id:
            no_neighbor_contigs += 1
            tot += 1
        else:
            tot += 1

    sample2no_n_contigs[sample] = f"{round(no_neighbor_contigs/tot*100, 2)}%"


mash_table = report.get_mash_table(mash_results, mash_detail, sample2scientific_name)
centrifuge_table = report.get_centrifuge_table(centrifuge_tables, sample2scientific_name)
rrna_table = report.get_rrna_summary_table(rrna_classification_file, sample2scientific_name)
multiqc_table = report.get_multiqc_table(assembly_multiqc=multiqc_assembly)
qualimap_table = report.qualimap_table(qualimap_reports, self_mapping=True)
checkm_table = report.checkm_table(checkm_table)


table_lowcoverage_contigs = quality_table(low_cov_fastas,
                                          sample2gc,
                                          sample2median_depth,
                                          sampls2cumulated_size,
                                          sampls2cumulated_size_filtered,
                                          sample2n_contigs,
                                          sample2no_n_contigs,
                                          sample2scientific_name,
                                          low_cov_detail=low_cov_detail,
                                          depth_cutoff=depth_cutoff)

table_resistance = resistance_table(resistance_reports)

report_str = f"""

.. raw:: html

    {SCRIPT}

    {STYLE}

=============================================================
Diag Pipeline - Resistance report
=============================================================

.. contents::
    :backlinks: none
    :depth: 2

Quality Control
---------------

MultiQC
*******

MultiQC aggregate results from bioinformatics analyses across many samples into a single report.
The analyses covered here include genome assembly with spades, evaluation of the sequencing
depth by mapping of the reads against the assembly and annotation with prokka.



.. raw:: html

    {multiqc_table}

Assembly overview
******************

.. raw:: html

    {table_lowcoverage_contigs}

Contamination check: mash (assembly)
************************************

Three best Mash hits (excluding phages).

.. raw:: html

    {mash_table}

Contamination check: Centrifuge
********************************

.. raw:: html

    {centrifuge_table}


Contamination check: 16S rRNA (contigs >500bp)
**********************************************

.. raw:: html

    {rrna_table}
    
Contamination check: checkM (contigs >500bp & depth >{depth_cutoff})
*********************************************************************

.. raw:: html

    {checkm_table}
    
Qualimap reports (self mapping)
********************************

.. raw:: html

    {qualimap_table}

Resistance (RGI/CARD)
----------------------


RGI Overview
*************

.. image:: {rgi_overview}
   :width: 60%

Detail
*******

.. raw:: html

    {table_resistance}


.. raw:: html

    <script>
    $(document).ready(function() {{
        $('#16S_table').DataTable( {{
            dom: 'Bfrtip',
            "pageLength": 10,
            "order": [[ 0, "desc" ]],
            "searching": true,
            "bLengthChange": false,
            "paging":   false,
            "info": false,
            'rowCallback': function(row, data, index){{
            var taxnonomy = data[4];
            var genus = data[2].split("_")[0];
            if (!(taxnonomy.includes(genus))){{
                            $(row).css('background-color', 'rgba(255, 0, 0, 0.7)');
                            }}
            }},
        }} );
    }} );
    </script>
"""
with open(output_file, "w") as fh:
    publish_file(
        source=io.StringIO(report_str),
        destination=fh,
        writer_name="html",
        settings_overrides={"stylesheet_path": ""},
    )
