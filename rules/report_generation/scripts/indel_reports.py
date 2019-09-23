
import pandas
import re
import numpy
import vcf
from Bio import SeqIO
# inputs
bed_file = snakemake.input["deletion_bed"]
gbk_file = snakemake.input["gbk_file"]
report_file = snakemake.output["report_file"]


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
    

def parse_gbk(gbk_file):
    record_dict = SeqIO.to_dict(SeqIO.parse(gbk_file, 'genbank'))
    return record_dict


def get_indel_orf(del_start, del_end, feature_list):
    '''
    Sarch overlap between bed ranges (region with depth of 0) and feature list.

    format of the bed file:
    contig_11	0	29	0
    contig_11	110	125	0
    contig_37	0	126	0
    '''
    match_type = False
    match_list = []
    for n, feature in enumerate(feature_list):
        if feature.type == 'source':
            continue
        if int(feature.location.start) < int(feature.location.start) and int(feature.location.end) < int(del_end):
            match_type = "overlap_begin"
        elif int(feature.location.start) >= int(del_start) and int(feature.location.end) <= int(del_end):
            match_type = "within"
        elif int(feature.location.start) > int(del_start) and int(feature.location.start) <= int(del_end) and int(feature.location.end) > int(del_end):
            match_type = "overlap_end"
        else:
            continue
        if match_type:
            match_list.append([feature, match_type])
            match_type = False
    return match_list


def bed2annotated_html_table(gbk_file, bed_file):
    records_dico = parse_gbk(gbk_file)

    bed_table = pandas.read_csv(bed_file, delimiter='\t', index_col=0, names=["record","start","end","depth"])
    gap_n = 0

    feature_table = []
    for i, row in bed_table.iterrows():
        gap_n += 1
        try:
            rec = records_dico[row.name]
        except KeyError:
            # contig removed based on length or low depth
            feature_table.append([row.name, "removed", gap_n, int(row[1]) - int(row[0]), "-", "-", "-", "-", "-"])
        else:

            match_list = get_indel_orf(row[0], row[1], rec.features)
            if len(match_list) > 0:
                for match in match_list:
                    if match[0].type == 'CDS':
                        if "gene" in match[0].qualifiers:
                            gene = match[0].qualifiers["gene"][0]
                        else:
                            gene = '-'
                    # of other kind of feature than CDS or gene
                    elif match[0].type != 'gene':
                        gene = '-'
                    else:
                        continue
                    feature_table.append([row.name,
                                          len(rec.seq),
                                          gap_n,
                                          int(row[1]) - int(row[0]),  # gap length
                                          match[1],
                                          match[0].location.start,
                                          match[0].location.end,
                                          match[0].type,
                                          gene])

            else:
                feature_table.append([row.name, len(rec.seq), gap_n, int(row[1]) - int(row[0]), "-", "-", "-", "-", "-"])

    header = ["record", "record length", "gap", "gap length", "localization", "feature_start", "feature_end", "type", "gene"]
    df = pandas.DataFrame(feature_table, columns=header)

    # cell content is truncated if colwidth not set to -1
    pandas.set_option('display.max_colwidth', -1)

    df_str = df.to_html(
        index=False,
        bold_rows=False,
        classes=["dataTable"],
        table_id="feature_table",
        escape=False,
        border=0)

    bed_table['contig'] = bed_table.index
    df_str2 = bed_table.to_html(
        index=False,
        bold_rows=False,
        classes=["dataTable"],
        table_id="gap_table",
        escape=False,
        border=0)

    return [df_str2.replace("\n", "\n" + 10 * " "),
            df_str.replace("\n", "\n" + 10 * " ")]

STYLE = """
    <link rel="stylesheet" type="text/css" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css"/>
    <link rel="stylesheet" type="text/css" href="https://cdn.datatables.net/1.10.19/css/dataTables.bootstrap.min.css"/>
    <style type="text/css">
    .docutils.container{width:100%}
    #gap_table {
        style="width:60%";
        word-break: normal;
    }
    #feature_table {
        style="width:60%";
        word-break: normal;
    }
    body{
      margin: 20px;
      padding: 20px;
    }
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
        $('#feature_table').DataTable( {
            dom: 'Bfrtip',
            "pageLength": 50,
            "searching": true,
            "bLengthChange": false,
            "paging":   true,
            "info": false,
            'rowCallback': function(row, data, index){
            if(data[16] == "-" ){
                             $(row).find('td:eq(16)').css('background-color', 'rgba(255, 0, 0, 0.6)');
                            }
            if(data[15] == "-" ){
                             $(row).find('td:eq(15)').css('background-color', 'rgba(255, 0, 0, 0.6)');
                            }
            if(data[14] == "-" ){
                             $(row).find('td:eq(14)').css('background-color', 'rgba(255, 0, 0, 0.6)');
                            }
            if(data[13] == "-" ){
                             $(row).find('td:eq(13)').css('background-color', 'rgba(255, 0, 0, 0.6)');
                            }
            if(data[11] == "YES"){
                             $(row).find('td:eq(11)').css('background-color', 'rgba(0, 255, 0, 0.6)');
                            }
            if(data[17] == "FAIL"){
                             $(row).find('td:eq(17)').css('background-color', 'rgba(255, 0, 0, 0.6)');
                            }
            if(data[17] == "PASS"){
                             $(row).find('td:eq(17)').css('background-color', 'rgba(0, 255, 0, 0.6)');
                            }
            if(data[16] == "True"){
                             $(row).find('td:eq(16)').css('background-color', 'rgba(255, 0, 0, 0.6)');
                            }
            },
            "columnDefs": [{
                            "targets": '_all',
                            "createdCell": function (td, cellData, rowData, row, col) {
                                $(td).css('padding-right', '20px')
                            }
                        }],
        } );
        $('#gap_table').DataTable( {
            dom: 'Bfrtip',
            "pageLength": 50,
            "searching": true,
            "bLengthChange": false,
            "paging":   true,
            "info": false,
            "columnDefs": [{
                            "targets": '_all',
                            "createdCell": function (td, cellData, rowData, row, col) {
                                $(td).css('padding-right', '20px')
                            }
                        }],
        } );

    } );
    </script>

    """


def write_report(output_file,
                 STYLE,
                 SCRIPT):

    import io
    from docutils.core import publish_file, publish_parts
    from docutils.parsers.rst import directives

    overview_table, feature_table = bed2annotated_html_table(gbk_file, bed_file)

    report_str = f"""

.. raw:: html

    {SCRIPT}

    {STYLE}

=============================================================
DELETIONs report
=============================================================


Overview Table
---------------

.. raw:: html


    {overview_table}

Feature Table
--------------

.. raw:: html


    {feature_table}

"""
    with open(output_file, "w") as fh:
        publish_file(
            source=io.StringIO(report_str),
            destination=fh,
            writer_name="html",
            settings_overrides={"stylesheet_path": ""},
        )

write_report(report_file,
             STYLE,
             SCRIPT)
