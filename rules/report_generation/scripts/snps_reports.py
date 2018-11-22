
import pandas
import re
import numpy

# inputs
vcf_file = snakemake.input["vcf_file"]
merged_vcf = snakemake.input["merged_vcf"]
gbk_file = snakemake.input["merged_vcf"]

# output
report_file = snakemake.output[0]


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


def parse_gbk(gbk_file):

    record_dict = SeqIO.to_dict(SeqIO.parse(gbk_file, 'genbank'))
    return record_dict


def get_neiboring_orf(position, feature_list):

    if int(position) > int(feature_list[-1].location.end):
        try:
            gene = feature_list[-1].qualifiers["gene"][0]
            locus_tag = feature_list[-1].qualifiers["locus_tag"][0]
        except KeyError:
            gene = '-'
            locus_tag = feature_list[-1].qualifiers["locus_tag"][0]
        return ['%s (%s)' % (locus_tag, gene), '-']
    first_feature = False
    for n, feature in enumerate(feature_list):
        if not first_feature:
            if int(position) < int(feature.location.start):

                try:
                    gene = feature.qualifiers["gene"][0]
                    locus_tag = feature.qualifiers["locus_tag"][0]
                except KeyError:
                    gene = '-'
                    locus_tag = feature.qualifiers["locus_tag"][0]
                return ["-", "%s (%s)" % (locus_tag, gene)]
        if feature.type != 'source':
            first_feature = True


        if int(position) > int(feature.location.end) and int(position) < int(feature_list[n + 1].location.start):
            try:
                gene1 = feature.qualifiers["gene"][0]
                locus_tag1 = feature.qualifiers["locus_tag"][0]
            except KeyError:
                gene1 = '-'
                locus_tag1 = feature.qualifiers["locus_tag"][0]
            try:
                gene2 = feature_list[n + 1].qualifiers["gene"][0]
                locus_tag2 = feature_list[n + 1].qualifiers["locus_tag"][0]
            except KeyError:
                gene2 = '-'
                locus_tag2 = feature_list[n + 1].qualifiers["locus_tag"][0]
            return ["%s (%s)" % (locus_tag1, gene1), "%s (%s)" % (locus_tag2, gene2)]
    return ["no CDS", "no CDS"]


def parse_vcf(vcf_file, gbk_file):
    from copy import copy
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import IUPAC
    from Bio.Seq import MutableSeq
    from Bio.Alphabet import generic_dna
    import numpy

    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    gbk_dico = parse_gbk(gbk_file)

    intergenic = 0
    synonymous = 0
    non_synonymous = 0
    indel = 0

    filter_head = ['%s (%s)' % (vcf_reader.filters[i].id, vcf_reader.filters[i].desc) for i in vcf_reader.filters ]
    header = ["contig", "contig_length", "position", "AC", "REF", "ALT", "MUTATION_location", "MUTATION_type", "ORF", "gene", "orf_before", "orf_after"] + filter_head
    table_rows = []
    for n, record in enumerate(vcf_reader):

        record_alt = copy(gbk_dico[record.CHROM])
        record_alt.seq = MutableSeq(str(record_alt.seq), generic_dna)
        # default to intergenic (modigied if feature match found)
        mut_location = 'Intergenic'
        mut_type = '-'
        orf_name = '-'
        gene = '-'
        for feature in record_alt.features:
            if int(record.POS) in feature and feature.type != "source":
                mut_location = feature.type
                orf_name = feature.qualifiers["locus_tag"][0]
                try:
                    gene = feature.qualifiers["gene"][0]
                except:
                    gene = '-'
                # print("SNP in feature--------------------------------------------:", feature.type, feature.location)
                if feature.type == 'CDS':

                    if len(record.ALT[0]) > 1:
                        mut_type = 'INDEL'
                        indel+=1
                        continue
                    else:
                        aa_seq_ref = str(feature.extract(record_alt.seq).translate())
                        # mutate reference sequence
                        record_alt.seq[int(record.POS)-1] = str(record.ALT[0])
                        # check if synonymous or not
                        aa_seq_alt = str(feature.extract(record_alt.seq).translate())

                        if str(aa_seq_ref) == str(aa_seq_alt):
                            mut_type = 'SYNONYMOUS'
                            synonymous += 1
                        else:
                            mut_type = 'NON SYNONYMOUS'
                            if len(record.FILTER) == 0:
                                non_synonymous+=1
                        #print (str(aa_seq_ref.translate()) == str(record_alt.translate()))
        if mut_location == 'Intergenic':
            orf_before, orf_after = get_neiboring_orf(int(record.POS), record_alt.features)
        else:
            orf_before, orf_after = ['-', '-']
        # print (filter_dico["PASS"].id)
        # print (filter_dico["PASS"].desc)
        contig = record.CHROM
        ac = record.INFO["AC"][0]
        position = record.POS
        ref = record.REF
        alt = record.ALT[0]
        filter_status = []

        # if any of the test failed, set PASS as failed
        if len(record.FILTER) != 0:
            record.FILTER.append('PASS')

        for filter_name in vcf_reader.filters:
            if filter_name in record.FILTER:
                if filter_name == 'PASS':
                    filter_status.append('NO')
                else:
                    filter_status.append('FAIL')
            else:
                if filter_name == 'PASS':
                    filter_status.append('YES')
                else:
                    filter_status.append('SUCCESS')

        table_rows.append("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (contig,
                                                                                  len(record_alt),
                                                                                  position,
                                                                                  ac,
                                                                                  ref,
                                                                                  alt,
                                                                                  mut_location,
                                                                                  mut_type,
                                                                                  orf_name,
                                                                                  gene,
                                                                                  orf_before,
                                                                                  orf_after,
                                                                                  '\t'.join(list(filter_status))))

    df = pandas.DataFrame(table_rows, columns=header)

    # cell content is truncated if colwidth not set to -1
    pandas.set_option('display.max_colwidth', -1)

    df_str = df.to_html(
        index=False,
        bold_rows=False,
        classes=["dataTable"],
        table_id="snps_table",
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
                 SCRIPT):

    import io
    from docutils.core import publish_file, publish_parts
    from docutils.parsers.rst import directives

    snp_table = parse_vcf(vcf_file, gbk_file)

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




Table
------

.. raw:: html

    <div class="alert alert-warning" role="alert">
      Genes with a hit <span class="label label-default">coverage &lt; 90%</span> are highlighted in <span class="label label-success">green</span> (if any) </br>
      Genes with an  <span class="label label-default">identity &lt; 90%</span> are highlighted in <span class="label label-danger">red</span> (if any)
    </div>

    {snp_table}

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
