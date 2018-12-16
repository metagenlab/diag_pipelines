
# from snakemake.utils import report

# with open(input[0]) as vcf:
#    n_calls = sum(1 for l in vcf if not l.startswith("#"))
# n_samples = list(read_naming.keys()
import pandas
from MN_tree import get_MN_tree, convert2cytoscapeJSON
from report import quality_table, plot_heatmap_snps, get_core_genome_size, get_reference_genome_size
import report
import io
from docutils.core import publish_file, publish_parts
from docutils.parsers.rst import directives

multiqc_assembly = snakemake.input["multiqc_assembly"]  # ok
contig_gc_depth_file_list = snakemake.input["contig_gc_depth_files"]
#snps_individual = snakemake.input["snps_individual"]

#### one per reference genome ###########################################
multiqc_mapping_list = snakemake.input["multiqc_mapping_list"]  # ok
snp_tables = snakemake.input["snp_tables"]  # ok
spanning_trees = snakemake.input["spanning_trees"]  # ok
reference_genomes = snakemake.input["reference_genomes"]  # ok
# multiple files of each reference genome
undetermined_snp_tables = snakemake.input["undetermined_positions"]  # ok
snps_reports = snakemake.input["snps_reports"]  # ok
indel_reports = snakemake.input["indel_reports"]
#snps_merged = snakemake.input["snps_merged"]
#########################################################################

ordered_samples = snakemake.params["samples"]

snp_detail_table = report.get_snp_detail_table(snps_reports, indel_reports)

# mlst_tree = snakemake.input["mlst_tree"]
mlst_tree = ""  #'/'.join(mlst_tree.split('/')[1:])

low_cov_fastas = snakemake.input["low_cov_fastas"]

# optional params (if cgMLST among the reference genomes)
core_genome_bed = snakemake.params["core_genome_bed"]

output_file = snakemake.output[0]

leaf2mlst = pandas.read_csv(snakemake.input["mlst"],
                           delimiter='\t',
                           names=["leaf","species","mlst","1","2","3","4","5","6","7"]).set_index("leaf").to_dict()["mlst"]


# get contig depth and GC
sample2gc = {}
sample2median_depth = {}
sampls2cumulated_size = {}
sample2n_contigs = {}
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
    sample2n_contigs[sample] = n_contigs

sample2scientific_name = snakemake.params["sample_table"].to_dict()["ScientificName"]

STYLE = """
    <link rel="stylesheet" type="text/css" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css"/>
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.19/css/dataTables.bootstrap.min.css"/>
    <style type="text/css">
    body{font-family:Helvetica,arial,sans-serif;font-size:14px;line-height:1.6;padding-bottom:10px;background-color:#fff;color:#333;margin:0}body>div .section::before{content:"";display:block;height:80px;margin:-80px 0 0}#quality-control::before{margin:0}.topic-title{font-size:18pt}body>div>.section{margin-left:255px;margin-bottom:3em}div.section{margin-right:20px}#contents>p{display:none}button,li p.first{display:inline-block}#contents{margin-top:80px;padding-left:0;width:235px;background-color:#f1f1f1;height:100%;position:fixed;overflow:auto}#contents ul{list-style-type:none}#contents ul>li{font-size:14pt}#contents ul>li a:hover{color:#151d26}button,h1.title{color:#fff;background-color:#151d26}#contents ul>li>ul>li{font-size:12pt}h1.title{margin-top:0;position:fixed;z-index:10;padding:20px;width:100%}code,table tr:nth-child(2n),tt{background-color:#f8f8f8}.one-col{min-width:310px;height:500px;margin:0 auto}.two-col-left{height:300px;width:49%;float:left}.two-col-right{height:300px;width:49%;float:right}button{margin:0 5px 0 0;padding:5px 25px;font-size:18px;line-height:1.8;appearance:none;box-shadow:none;border-radius:3px;border:none}button:focus{outline:0}button:hover{background-color:#4183C4}button:active{background-color:#27496d}.legend-rect{width:20px;height:20px;margin-right:8px;margin-left:20px;float:left;-webkit-border-radius:2px;border-radius:2px}a{color:#4183C4;text-decoration:none}a.absent{color:#c00}a.anchor{padding-left:30px;margin-left:-30px;cursor:pointer;position:absolute;top:0;left:0;bottom:0}dl,dl dt,dl dt:first-child,hr,table,table tr{padding:0}table tr td,table tr th{border:1px solid #ccc;text-align:left;padding:6px 13px}h1,h2,h3,h4,h5,h6{margin:20px 0 10px;padding:0;font-weight:700;-webkit-font-smoothing:antialiased;cursor:text;position:relative}h1:hover a.anchor,h2:hover a.anchor,h3:hover a.anchor,h4:hover a.anchor,h5:hover a.anchor,h6:hover a.anchor{text-decoration:none}h1 code,h1 tt,h2 code,h2 tt,h3 code,h3 tt,h4 code,h4 tt,h5 code,h5 tt,h6 code,h6 tt{font-size:inherit}h1{font-size:28px;color:#151d26;border-bottom:1px solid #ccc}h2{font-size:24px;color:#000}h3{font-size:18px}h4{font-size:16px}dl dt,h5,h6{font-size:14px}h6{color:#777}blockquote,dl,li,ol,p,pre,table,ul{margin:15px 0}hr{background:url(http://tinyurl.com/bq5kskr) repeat-x;border:0;color:#ccc;height:4px}a:first-child h1,a:first-child h2,a:first-child h3,a:first-child h4,a:first-child h5,a:first-child h6{margin-top:0;padding-top:0}h1 p,h2 p,h3 p,h4 p,h5 p,h6 p{margin-top:0}dl dt{font-weight:700;font-style:italic;margin:15px 0 5px}blockquote>:first-child,dl dd>:first-child,dl dt>:first-child,table tr td :first-child,table tr th :first-child{margin-top:0}blockquote>:last-child,dl dd>:last-child,dl dt>:last-child{margin-bottom:0}dl dd{margin:0 0 15px;padding:0 15px}blockquote{border-left:4px solid #ddd;padding:0 15px;color:#777}table{border-spacing:0;border-collapse:collapse}table tr{border-top:1px solid #ccc;background-color:#fff;margin:0}table tr th{font-weight:700;margin:0}table tr td{margin:0}table tr td :last-child,table tr th :last-child{margin-bottom:0}img{max-width:100%}span.frame{display:block;overflow:hidden}span.frame>span{border:1px solid #ddd;display:block;float:left;overflow:hidden;margin:13px 0 0;padding:7px;width:auto}span.frame span img{display:block;float:left}span.frame span span{clear:both;color:#333;display:block;padding:5px 0 0}span.align-center{display:block;overflow:hidden;clear:both}span.align-center>span{display:block;overflow:hidden;margin:13px auto 0;text-align:center}span.align-center span img{margin:0 auto;text-align:center}span.align-right{display:block;overflow:hidden;clear:both}span.align-right>span{display:block;overflow:hidden;margin:13px 0 0;text-align:right}span.align-right span img{margin:0;text-align:right}span.float-left{display:block;margin-right:13px;overflow:hidden;float:left}span.float-left span{margin:13px 0 0}span.float-right{display:block;margin-left:13px;overflow:hidden;float:right}span.float-right>span{display:block;overflow:hidden;margin:13px auto 0;text-align:right}code,tt{margin:0 2px;padding:0 5px;white-space:nowrap;border:1px solid #eaeaea;border-radius:3px}pre code{margin:0;padding:0;white-space:pre;background:0 0}.highlight pre,pre{background-color:#f8f8f8;border:1px solid #ccc;font-size:13px;line-height:19px;overflow:auto;padding:6px 10px;border-radius:3px}pre code,pre tt{background-color:transparent;border:none}div#metadata{text-align:right}h1{line-height:1.6}.simple{padding-left:20px}.docutils.container{width:100%}
    .pull-left{
    .dataTables_filter {
       float: left !important;
    }
    #cy {
        width: 100%;
        height: 100%;
        position: absolute;
        left: 0;
        top: 0;
        z-index: 999;
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
        $('#cov_table').DataTable( {
            dom: 'Bfrtip',
            "pageLength": 20,
            "searching": true,
            "bLengthChange": false,
            "paging":   true,
            "info": false
        } );
    } );

    </script>

    """

if core_genome_bed:
    core_size = get_core_genome_size(core_genome_bed)
    for i in reference_genomes:
        if "cgMLST" in i:
            ref_size = get_reference_genome_size(i)
    fraction_core = round(float(core_size) / float(ref_size) * 100, 2)
    core_str = """
    - Size of the reference genome: %s
    - Size of the core genome: %s (%s %% of the reference)
    """ % (ref_size,
           core_size,
           fraction_core)
else:
    core_size = False
    core_str = ""

multiqc_table = report.get_multiqc_table(multiqc_assembly,
                                         multiqc_mapping_list)

table_lowcoverage_contigs = quality_table(low_cov_fastas,
                                          sample2gc,
                                          sample2median_depth,
                                          sampls2cumulated_size,
                                          sample2n_contigs,
                                          sample2scientific_name,
                                          undetermined_snps_files=undetermined_snp_tables,
                                          core_genome_size=core_size)

snp_heatmap_str = ""
for n, snp_table in enumerate(snp_tables):
    snp_heatmap_str += '''

%s
******************************************************

.. raw:: html

    %s


    ''' % (snp_table.split("/")[2],
           plot_heatmap_snps(snp_table, id=snp_table.split("/")[2] ))

spanning_tree_str = ""
for n, tree in enumerate(spanning_trees):
    tree_path = '/'.join(tree.split('/')[1:])
    spanning_tree_str += """

%s
***************************************************************
.. figure:: %s
   :alt: %s
   :width: 60%%
    """ % (tree.split("/")[3], tree_path, tree.split("/")[3])

report_str = f"""

.. raw:: html

    {SCRIPT}

    {STYLE}

=============================================================
Diag Pipeline - Epidemiology
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
If reads were mapped against multiple reference genomes, on report per reference is generated.


.. raw:: html

    {multiqc_table}

Overview quality
*****************

.. raw:: html

    {table_lowcoverage_contigs}

Typing
------

MS tree(s)
-----------

{spanning_tree_str}


SNP table(s)
-------------

{snp_heatmap_str}


Variant details
***************

Detailed list of variants identified compared to one or multiple reference genomes.

.. raw:: html

    {snp_detail_table}

"""

"""

Phylogeny + MLST
****************

.. figure:: {mlst_tree}
   :alt: MST tree
   :width: 40%

   MLST as determined by T. Seemann mlst_.

"""

with open(output_file, "w") as fh:
    publish_file(
        source=io.StringIO(report_str),
        destination=fh,
        writer_name="html",
        settings_overrides={"stylesheet_path": ""},
    )



#net = ''#convert2cytoscapeJSON(get_MN_tree(snp_table), leaf2mlst)
