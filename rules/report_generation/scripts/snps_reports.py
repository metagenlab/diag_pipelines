
import pandas
import re
import numpy
import vcf
from Bio import SeqIO

# inputs
vcf_file = snakemake.input["vcf_file"]
merged_vcf = snakemake.input["merged_vcf"]
gbk_file = snakemake.input["gbk_file"]
reference = snakemake.params["reference"]
if "_assembled_genome" in reference:
    reference = re.sub("_assembled_genome", "", reference)
    assembleg_genome = True

# output
report_file = snakemake.output["html_file"]

merged_vcf_records = [i for i in vcf.Reader(open(merged_vcf, 'r'))]

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
    print("gb", gbk_file)
    record_dict = SeqIO.to_dict(SeqIO.parse(gbk_file, 'genbank'))
    print("record_dict", record_dict)

    return record_dict


def get_neiboring_orf(position, feature_list):
    '''
    Identify neiboring feature of variant position.
    Input: SeqRecord and position (integer)
    Output: List of two strings with the closest features located before and after the input position
    (either locus_tags/gene names).
    If no feature found (eg. in the case where the locate in close to the begining/end of the record), return "no CDS".
    '''
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


def check_reference_mapping_GT(vcf_record_list,
                               contig,
                               position,
                               sampe_reference):
    '''
    Quality check. Check if a variant position is also a variant in the mapping
    the mapping of the reference reads versus the reference assembly.
    Only necessary when dealing with assembled references.
    Input: merged VCF file with all variant positions (all_samples_snp.vcf).
    Return False if a variant is found at the given position.
    '''
    print("vcf check--------------")
    print(contig, position, sampe_reference)
    for record in vcf_record_list:
        if record.CHROM == contig and record.POS == position:
            # search for reference sample
            for sample in record.samples:
                if sample.sample == sampe_reference:
                    # check if ALT position more supported than REF
                    if len(record.FILTER) == 0:
                        PASS = "PASS"
                    else:
                        PASS = "FAIL"
                    print(sample.sample)
                    print(type(sample['GT']))
                    if sample['GT'] not in ['.', "0"]:
                        return [True, PASS]
                    else:
                        return [False, PASS]
    return [None, None]

def parse_vcf(vcf_file, gbk_file):
    '''
    Given a vcf input file and the gbk of the reference genome, return an html table of identified variants.
    '''

    from copy import copy
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import IUPAC
    from Bio.Seq import MutableSeq
    from Bio.Alphabet import generic_dna
    import numpy

    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    gbk_dico = parse_gbk(gbk_file)
    print(gbk_dico.keys())
    intergenic = 0
    synonymous = 0
    non_synonymous = 0
    indel = 0

    filter_head = ['%s' % (vcf_reader.filters[i].id) for i in vcf_reader.filters]
    header = ["contig", "length", "position", "REF", "ALT", "location", "type", "ORF", "gene", "orf_before", "orf_after"]
    header += filter_head
    if 'assembled' in vcf_file:
        header.append("InRef")
        header.append("Fail Others")

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
                if feature.type == 'mobile_element':
                    orf_name = feature.qualifiers["mobile_element_type"][0]
                elif feature.type == 'CDS':
                    orf_name = feature.qualifiers["locus_tag"][0]
                else:
                    orf_name = "Unknown locus for feature: %s" % feature.type
                try:
                    gene = feature.qualifiers["gene"][0]
                except:
                    gene = '-'
                # print("SNP in feature--------------------------------------------:", feature.type, feature.location)
                if feature.type == 'CDS':

                    if len(record.ALT[0]) > 1:
                        mut_type = 'INDEL'
                        indel += 1
                        continue
                    else:
                        aa_seq_ref = str(feature.extract(record_alt.seq).translate())
                        # mutate reference sequence
                        record_alt.seq[int(record.POS) - 1] = str(record.ALT[0])
                        # check if synonymous or not
                        aa_seq_alt = str(feature.extract(record_alt.seq).translate())

                        if str(aa_seq_ref) == str(aa_seq_alt):
                            mut_type = 'S'
                            synonymous += 1
                        else:
                            mut_type = 'NS'
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
        print("depth", record.samples[0]['AD'], record.samples[0]['DP'])

        if record.samples[0]['GT'] in ['.']:
            # print(type(ac))
            # if ac == 0:
            continue
        position = record.POS
        ref = "%s (%s/%s)" % (record.REF,
                              record.samples[0]['AD'][0],
                              record.samples[0]['DP'])
        if len(record.ALT[0]) == 1:
            alt = "%s (%s/%s)" % (record.ALT[0],
                                  record.samples[0]['AD'][1],
                                  record.samples[0]['DP'])
        else:
            alt = "%sbp (%s/%s)" % (len(record.ALT[0]),
                                    record.samples[0]['AD'][1],
                                    record.samples[0]['DP'])
        filter_status = []

        # if any of the test failed, set PASS as failed
        if len(record.FILTER) != 0:
            record.FILTER.append('PASS')

        for filter_name in vcf_reader.filters:
            if filter_name in record.FILTER:
                if filter_name == 'PASS':
                    filter_status.append('NO')
                else:
                    filter_status.append('-')
            else:
                if filter_name == 'PASS':
                    filter_status.append('YES')
                else:
                    filter_status.append('+')
        row = [contig,
               len(record_alt),
               position,
               ref,
               alt,
               mut_location,
               mut_type,
               orf_name,
               gene,
               orf_before,
               orf_after]

        row += list(filter_status)

        if 'assembled' in vcf_file:
            GT, PASS = check_reference_mapping_GT(merged_vcf_records,
                                            contig,
                                            position,
                                            reference)
            print("################# GT #################")
            row.append(GT)
            row.append(PASS)

        table_rows.append(row)

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
    .docutils.container{width:100%}
    #snps_table {
        style="width:100%";
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
        $('#snps_table').DataTable( {
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
SNPS report
=============================================================

Mapping on raw assembly, but filter based on gbk file annotted with prokka after filtering of small and spurious contigs (SNPS identified on spurious contigs, e.g low coverage contigs, are not reported in this table).

Table
------

===========  ===================================================================================================
Header name  Description
===========  ===================================================================================================
Contig       Contig name
Length       Length of the contig
Position     Position in the contig
REF          Reference
ALT          Alternative variant
Location     Localization within or between ORFs
Type         Either Indel, Synonymous mutation OR non synonymous mutation
ORF          ORF accession
Gene         Gene name
orf_before   Closest ORF locate before the variant (intergenic only)
orf_after    Closest ORF located after the variant (intergenic only)
PASS         Passed all quality filters?
Lowqual      Result of the quality filter (quality-score-threshold: 18)
Lowcov       Result of the depth filter (fail if depth < 10)
freqalt      Filter based on the minimum frequency of the alternative variant (fail if REF & REF freq < 75%)
freqref      Filter based on the minimum frequency of the reference variant (fail if ALT & ALT freq < 75%)
inref        Is this variant also present when mapping reads used to assemble the reference assembly? (optional)
fail others  Is this variant also found in other sample but failed quality filters?
===========  ===================================================================================================


.. raw:: html


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
