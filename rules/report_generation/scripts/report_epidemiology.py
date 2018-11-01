
# from snakemake.utils import report

# with open(input[0]) as vcf:
#    n_calls = sum(1 for l in vcf if not l.startswith("#"))
# n_samples = list(read_naming.keys()
import pandas
from MN_tree import get_MN_tree, convert2cytoscapeJSON
from report import coverage_table, virulence_table, resistance_table, plot_heatmap_snps, get_core_genome_size, get_reference_genome_size

multiqc_report = snakemake.input["multiqc_report"]
snp_table = snakemake.input["snp_table"]

ordered_samples = snakemake.params["samples"]
spanning_tree_core = snakemake.input["spanning_tree_core"]
spanning_tree_core = '/'.join(spanning_tree_core.split('/')[1:])

mlst_tree = snakemake.input["mlst_tree"]
mlst_tree = '/'.join(mlst_tree.split('/')[1:])

low_cov_fastas = snakemake.input["low_cov_fastas"]

reference_genome = snakemake.input["reference_genome"]
core_genome_bed = snakemake.input["core_genome_bed"]

output_file = snakemake.output[0]

leaf2mlst= pandas.read_csv(snakemake.input["mlst"],
                           delimiter='\t',
                           names=["leaf","species","mlst","1","2","3","4","5","6","7"]).set_index("leaf").to_dict()["mlst"]


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

    document.addEventListener('DOMContentLoaded', function(){
    
        var cy = window.cy = cytoscape({
            container: document.getElementById('cy'),
    
            autounselectify: true,
            
            boxSelectionEnabled: false,
    
            layout: {
                name: 'cose',
                idealEdgeLength: function(edge){ return 1/Math.log(edge.data('strength')); }, // edgeLength
                padding: 30,
                maxSimulationTime: 6000,
                randomize: false,
                nodeSpacing: 100,
                animate: true
    
            },
    
            style: [
                {
                    selector: 'node',
                    css: {
                        'background-color': 'data(color)',
                        'shape': 'roundrectangle',
                        'width': function(node){ return 3*node.data('label').length; },
                        'height': 7,
                        'content': 'data(label)',
                        'font-size': 5,
                        'text-valign': 'center'
                    }
                },
    
                {
                    selector: 'edge',
                    css: {
                        'line-color': 'data(color)',
                        'label': 'data(strength)',
                        'width': 'data(width)',
                        'font-size': 4
                    }
                }
            ],
    
            elements:  %s
            });
    
    });
    </script>

    """



def write_report(output_file,
                 STYLE,
                 SCRIPT,
                 spanning_tree_core,
                 low_cov_fasta,
                 mlst_tree,
                 snp_table,
                 core_genome_size,
                 reference_genome_size):
    import io
    from docutils.core import publish_file, publish_parts
    from docutils.parsers.rst import directives

    multiqc_link = '<a href="%s">MiltiQC</a>' % '/'.join(multiqc_report.split('/')[1:])
    table_lowcoverage_contigs = coverage_table(low_cov_fasta)

    snp_heatmap = plot_heatmap_snps(snp_table)
    fraction_core = round(float(core_genome_size)/float(reference_genome_size)*100, 2)

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


.. raw:: html

    {multiqc_link}
    
Low coverage contigs
********************

.. raw:: html

    {table_lowcoverage_contigs}
    
Typing
------

MLST
*****

The *S. aureus* MLST scheme is based on the sequence of the following seven house-keeping genes:
    
1. arcC (Carbamate kinase)
2. aroE (Shikimate dehydrogenase)
3. glpF (Glycerol kinase)
4. gmk (Guanylate kinase)
5. pta (Phosphate acetyltransferase)
6. tpi (Triosephosphate isomerase)
7. yqi (Acetyle coenzyme A acetyltransferase)
              
The MLST was determined using the mlst_ software based on PubMLST_ typing schemes.
    
.. _PubMLST: https://pubmlst.org/
.. _mlst: https://github.com/tseemann/mlst

Phylogeny + MLST
****************

.. figure:: {mlst_tree} 
   :alt: MST tree
   :figwidth: 80%

   This is the caption of the figure (a simple paragraph).

MS tree (R)
*********************

- Size of the reference genome: {reference_genome_size}
- Size of the core genome: {core_genome_size} ({fraction_core} % of the reference) 


.. figure:: {spanning_tree_core} 
   :alt: MST tree
   :figwidth: 80%

   This is the caption of the figure (a simple paragraph).
   
MS tree (js) 
***********************

.. raw:: html

    <div id="cy" style="width:80%;height:700px; position: relative; border: 2px solid #212523"></div>

SNP table
***********

.. raw:: html

    {snp_heatmap}
    
"""
    with open(output_file, "w") as fh:
        publish_file(
            source=io.StringIO(report_str),
            destination=fh,
            writer_name="html",
            settings_overrides={"stylesheet_path": ""},
        )



net = convert2cytoscapeJSON(get_MN_tree(snp_table), leaf2mlst)

write_report(output_file,
             STYLE,
             SCRIPT % net,
             spanning_tree_core,
             low_cov_fastas,
             mlst_tree,
             snp_table,
             get_core_genome_size(core_genome_bed),
             get_reference_genome_size(reference_genome))