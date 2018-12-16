
import pandas
from Bio import SeqIO
import report
from Bio.SeqUtils import GC
import io
from docutils.core import publish_file, publish_parts
import re

link = snakemake.input[0]
output_file = snakemake.output[0]
expected_genus = snakemake.params["expected_genus"]

def get_mash_table(link):
    header = ["score",
              "sketch",
               "n",
              "e-value",
              "accession",
              "description"]

    table = pandas.read_csv(link, delimiter="\t", names=header)
    table = table.sort_values(table.columns[0], ascending=False)
    rows = []
    for n, row in table.iterrows():
        score = row["score"]
        sketch = row["sketch"]
        evalue = row["e-value"]
        accession = re.findall("(GCF_[0-9]*\.[0-9]+)_.*", row["accession"])[0]
        description = re.sub("\[[0-9]+ seqs\] ", "", row["description"])


        data = [score,
                sketch,
                evalue,
                '<a href="https://www.ncbi.nlm.nih.gov/assembly/%s/">%s</a>' % (accession, accession),
                ' '.join(description.split(" ")[1:])[0:80]]
        rows.append(data)
    header = ["score",
              "sketch",
              "e-value",
              "accession",
              "description"]
    df = pandas.DataFrame(rows, columns=header)

    # cell content is truncated if colwidth not set to -1
    pandas.set_option('display.max_colwidth', -1)

    df_str = df.to_html(
        index=False,
        bold_rows=False,
        classes=["dataTable"],
        table_id="mash_table",
        escape=False,
        border=0)

    return df_str.replace("\n", "\n" + 10 * " ")


STYLE = """
    <link rel="stylesheet" type="text/css" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css"/>
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.19/css/dataTables.bootstrap.min.css"/>
    <style type="text/css">
    body{font-family:Helvetica,arial,sans-serif;font-size:14px;line-height:1.6;padding-bottom:10px;background-color:#fff;color:#333;margin:0}body>div .section::before{content:"";display:block;height:80px;margin:-80px 0 0}#mash-results::before{margin:0}.topic-title{font-size:18pt}body>div>.section{margin-left:255px;margin-bottom:3em}div.section{margin-right:20px}#contents>p{display:none}button,li p.first{display:inline-block}#contents{margin-top:80px;padding-left:0;width:235px;background-color:#f1f1f1;height:100%;position:fixed;overflow:auto}#contents ul{list-style-type:none}#contents ul>li{font-size:14pt}#contents ul>li a:hover{color:#151d26}button,h1.title{color:#fff;background-color:#151d26}#contents ul>li>ul>li{font-size:12pt}h1.title{margin-top:0;position:fixed;z-index:10;padding:20px;width:100%}code,table tr:nth-child(2n),tt{background-color:#f8f8f8}.one-col{min-width:310px;height:500px;margin:0 auto}.two-col-left{height:300px;width:49%;float:left}.two-col-right{height:300px;width:49%;float:right}button{margin:0 5px 0 0;padding:5px 25px;font-size:18px;line-height:1.8;appearance:none;box-shadow:none;border-radius:3px;border:none}button:focus{outline:0}button:hover{background-color:#4183C4}button:active{background-color:#27496d}.legend-rect{width:20px;height:20px;margin-right:8px;margin-left:20px;float:left;-webkit-border-radius:2px;border-radius:2px}a{color:#4183C4;text-decoration:none}a.absent{color:#c00}a.anchor{padding-left:30px;margin-left:-30px;cursor:pointer;position:absolute;top:0;left:0;bottom:0}dl,dl dt,dl dt:first-child,hr,table,table tr{padding:0}table tr td,table tr th{border:1px solid #ccc;text-align:left;padding:6px 13px}h1,h2,h3,h4,h5,h6{margin:20px 0 10px;padding:0;font-weight:700;-webkit-font-smoothing:antialiased;cursor:text;position:relative}h1:hover a.anchor,h2:hover a.anchor,h3:hover a.anchor,h4:hover a.anchor,h5:hover a.anchor,h6:hover a.anchor{text-decoration:none}h1 code,h1 tt,h2 code,h2 tt,h3 code,h3 tt,h4 code,h4 tt,h5 code,h5 tt,h6 code,h6 tt{font-size:inherit}h1{font-size:28px;color:#151d26;border-bottom:1px solid #ccc}h2{font-size:24px;color:#000}h3{font-size:18px}h4{font-size:16px}dl dt,h5,h6{font-size:14px}h6{color:#777}blockquote,dl,li,ol,p,pre,table,ul{margin:15px 0}hr{background:url(http://tinyurl.com/bq5kskr) repeat-x;border:0;color:#ccc;height:4px}a:first-child h1,a:first-child h2,a:first-child h3,a:first-child h4,a:first-child h5,a:first-child h6{margin-top:0;padding-top:0}h1 p,h2 p,h3 p,h4 p,h5 p,h6 p{margin-top:0}dl dt{font-weight:700;font-style:italic;margin:15px 0 5px}blockquote>:first-child,dl dd>:first-child,dl dt>:first-child,table tr td :first-child,table tr th :first-child{margin-top:0}blockquote>:last-child,dl dd>:last-child,dl dt>:last-child{margin-bottom:0}dl dd{margin:0 0 15px;padding:0 15px}blockquote{border-left:4px solid #ddd;padding:0 15px;color:#777}table{border-spacing:0;border-collapse:collapse}table tr{border-top:1px solid #ccc;background-color:#fff;margin:0}table tr th{font-weight:700;margin:0}table tr td{margin:0}table tr td :last-child,table tr th :last-child{margin-bottom:0}img{max-width:100%}span.frame{display:block;overflow:hidden}span.frame>span{border:1px solid #ddd;display:block;float:left;overflow:hidden;margin:13px 0 0;padding:7px;width:auto}span.frame span img{display:block;float:left}span.frame span span{clear:both;color:#333;display:block;padding:5px 0 0}span.align-center{display:block;overflow:hidden;clear:both}span.align-center>span{display:block;overflow:hidden;margin:13px auto 0;text-align:center}span.align-center span img{margin:0 auto;text-align:center}span.align-right{display:block;overflow:hidden;clear:both}span.align-right>span{display:block;overflow:hidden;margin:13px 0 0;text-align:right}span.align-right span img{margin:0;text-align:right}span.float-left{display:block;margin-right:13px;overflow:hidden;float:left}span.float-left span{margin:13px 0 0}span.float-right{display:block;margin-left:13px;overflow:hidden;float:right}span.float-right>span{display:block;overflow:hidden;margin:13px auto 0;text-align:right}code,tt{margin:0 2px;padding:0 5px;white-space:nowrap;border:1px solid #eaeaea;border-radius:3px}pre code{margin:0;padding:0;white-space:pre;background:0 0}.highlight pre,pre{background-color:#f8f8f8;border:1px solid #ccc;font-size:13px;line-height:19px;overflow:auto;padding:6px 10px;border-radius:3px}pre code,pre tt{background-color:transparent;border:none}div#metadata{text-align:right}h1{line-height:1.6}.simple{padding-left:20px}.docutils.container{width:100%}
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

    <script src="https://cdn.datatables.net/buttons/1.5.2/js/dataTables.buttons.min.js"></script>
    <script src="https://cdn.datatables.net/buttons/1.5.2/js/buttons.flash.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/jszip/3.1.3/jszip.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.1.36/pdfmake.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.1.36/vfs_fonts.js"></script>
    <script src="https://cdn.datatables.net/buttons/1.5.2/js/buttons.html5.min.js"></script>



    <script>
    $(document).ready(function() {
        $('#mash_table').DataTable( {
            dom: 'Bfrtip',
            buttons: [
            'copy', 'csv', 'excel', 'pdf'
            ],
            "pageLength": 20,
            "order": [[ 0, "desc" ]],
            "searching": true,
            "bLengthChange": false,
            "paging":   true,
            "info": false,
            'rowCallback': function(row, data, index){
            var string = data[4];
            if(string.includes("phage")){
                             $(row).css('background-color', 'rgba(200, 222, 222, 0.8)');
                            }
            if(string.includes("plasmid")){
                             $(row).css('background-color', 'rgba(200, 100, 222, 0.7)');
                            }
            if (!(string.includes("%s"))){
                             $(row).css('background-color', 'rgba(255, 0, 0, 0.7)');
                            }
            },
        } );
    } );
    </script>

    """ % expected_genus




contig_table = get_mash_table(link)

report_str = f"""

.. raw:: html

    {SCRIPT}

    {STYLE}

=============================================================
Diag Pipeline - Mash results
=============================================================

.. contents::
    :backlinks: none
    :depth: 2

Mash results
---------------------------

.. raw:: html

    {contig_table}

"""
with open(output_file, "w") as fh:
    publish_file(
        source=io.StringIO(report_str),
        destination=fh,
        writer_name="html",
        settings_overrides={"stylesheet_path": ""},
    )
