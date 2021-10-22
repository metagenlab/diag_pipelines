
#!/usr/bin/env python

import pandas
import re
import vcf
from Bio import SeqIO
import io
from docutils.core import publish_file
import codecs

# inputs
vcf_file = snakemake.input["vcf_file"]
merged_vcf = snakemake.input["merged_vcf"]
gbk_file = snakemake.input["gbk_file"]
reference = snakemake.params["reference"]

# rename reference if assembled genome
if "_assembled_genome" in reference:
    reference = re.sub("_assembled_genome", "", reference)

# output
report_file = snakemake.output["html_file"]

# parse vcf
merged_vcf_records = [i for i in vcf.Reader(codecs.open(merged_vcf, 'r', 'latin-1'))]


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
    with open(gbk_file, "r") as f:
        record_dict = SeqIO.to_dict(SeqIO.parse(f, 'genbank'))
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
        if 'gene' in feature_list[-1].qualifiers:
            gene = feature_list[-1].qualifiers["gene"][0]
            
        else:
            gene = '-'
        if 'locus_tag' in feature_list[-1].qualifiers:
            locus_tag = feature_list[-1].qualifiers["locus_tag"][0]
        else:
            locus_tag = '-'
        if 'mobile_element_type' in feature_list[-1].qualifiers:
            gene = feature_list[-1].qualifiers["mobile_element_type"]
        else:
            gene = '-'
        return ['%s (%s)' % (locus_tag, gene), '-']

    for n, feature in enumerate(feature_list):
        if feature.type in ['source', 'regulatory']:
            continue
        if int(position) < int(feature.location.start):
            if 'gene' in feature.qualifiers:
                gene = feature.qualifiers["gene"][0]
                locus_tag = feature.qualifiers["locus_tag"][0]
            else:
                gene = '-' 
            if 'locus_tag' in feature.qualifiers:
                locus_tag = feature.qualifiers["locus_tag"][0]
            if 'mobile_element_type' in feature.qualifiers:
                gene = feature.qualifiers["mobile_element_type"]
            else:
                gene = '-'
            return ["-", "%s (%s)" % (locus_tag, gene)]

        if int(position) > int(feature.location.end) and int(position) < int(feature_list[n + 1].location.start):
            if 'gene' in feature.qualifiers:
                gene1 = feature.qualifiers["gene"][0]
            else:
                gene1 = '-'
            if 'locus_tag' in feature.qualifiers:
                locus_tag1 = feature.qualifiers["locus_tag"][0]
            else:
                locus_tag1 = '-'
            if 'mobile_element_type' in feature.qualifiers:
                gene1 = feature.qualifiers["mobile_element_type"]
            else:
                gene1 = '-'
            if 'gene' in feature_list[n + 1].qualifiers:
                gene2 = feature_list[n + 1].qualifiers["gene"][0]
            else:
                gene2 = '-'
            if 'locus_tag' in feature_list[n + 1].qualifiers:
                locus_tag2 = feature_list[n + 1].qualifiers["locus_tag"][0]
            else:
                locus_tag2 = '-'
            if 'mobile_element_type' in feature_list[n + 1].qualifiers:
                gene2 = feature_list[n + 1].qualifiers["mobile_element_type"]
            else:
                gene2 = '-'
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
                    if sample['GT'] not in ['.', "0"]:
                        return [True, PASS]
                    else:
                        return [False, PASS]
    return [None, None]


def extract_mutation(seq1,seq2):
    for n, (aa1,aa2) in enumerate(zip(list(seq1), list(seq2))):
        if aa1 != aa2:
            return f'{aa1}{n+1}{aa2}'


def search_mutated_feature(vcf_record, gbk_dico):
    '''
    - Search if mutation is located within a coding sequence
    - determine if mutation is synonymous or not using a MutableSeq record (copy of the original record with mutation)
    '''
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from copy import copy
    from Bio.Alphabet import IUPAC
    from Bio.Seq import MutableSeq
    from Bio.Alphabet import generic_dna

    # create
    record_alt = copy(gbk_dico[vcf_record.CHROM])
    record_alt.seq = MutableSeq(str(record_alt.seq), generic_dna)

    results = {"mut_location": "Intergenic",
               "mut_type": '-',
               "orf_name": '-',
               "gene": '-'}

    for feature in record_alt.features:
        if int(vcf_record.POS) in feature and feature.type != "source":
            results["mut_location"] = feature.type
            if feature.type == 'mobile_element':
                results["orf_name"] = feature.qualifiers["mobile_element_type"][0]
            elif feature.type == 'CDS':
                results["orf_name"] = feature.qualifiers["locus_tag"][0]
            else:
                results["orf_name"] = "Unknown locus for feature: %s" % feature.type
            try:
                results["gene"] = feature.qualifiers["gene"][0]
            except KeyError:
                results["gene"] = '-'
            if feature.type == 'CDS':

                if len(vcf_record.ALT[0]) > 1:
                    results["mut_type"] = 'INDEL'
                    continue
                else:
                    aa_seq_ref = str(feature.extract(record_alt.seq).translate())
                    # mutate reference sequence
                    if vcf_record.ALT[0] == '*':
                    # frameshift
                        results["mut_type"] = 'F'
                    else:
                        record_alt.seq[int(vcf_record.POS) - 1] = str(vcf_record.ALT[0])

                        # check if synonymous or not
                        aa_seq_alt = str(feature.extract(record_alt.seq).translate())
                        if str(aa_seq_ref) == str(aa_seq_alt):
                            results["mut_type"] = 'S'
                        else:
                            results["mut_type"] = extract_mutation(aa_seq_ref,
                                                                   aa_seq_alt)
                    
            return results
    # if no match, return empty results
    return results


def parse_vcf(vcf_file, gbk_file):
    print(gbk_file)
    '''
    Given a vcf input file and the gbk of the reference genome, return an html table of identified variants.
    '''
    vcf_reader = vcf.Reader(codecs.open(vcf_file, 'r', 'latin-1'))
    gbk_dico = parse_gbk(gbk_file)

    filter_head = ['%s' % (vcf_reader.filters[i].id) for i in vcf_reader.filters]
    header = ["contig", "length", "position", "REF", "ALT", "location", "type", "ORF", "gene", "orf_before", "orf_after"]
    header += filter_head
    if 'assembled' in vcf_file:
        header.append("InRef")
        header.append("Fail Others")

    table_rows = []
    snp_count = 0
    for n, vcf_record in enumerate(vcf_reader):
        try:
            contig = gbk_dico[vcf_record.CHROM]
        except KeyError:
            print("Missing contig", vcf_record.CHROM)
            continue
        variant_feature = search_mutated_feature(vcf_record, gbk_dico)

        if variant_feature["mut_location"] == 'Intergenic':
            orf_before, orf_after = get_neiboring_orf(int(vcf_record.POS), contig.features)
        else:
            orf_before, orf_after = ['-', '-']

        contig_name = vcf_record.CHROM

        # skip ppositions with genomtype identical to REF
        if vcf_record.samples[0]['GT'] in ['.', '0']:
            continue
        snp_count += 1
        position = vcf_record.POS

        #  REF and ALT with respective depth in parenthesis
        ref = "%s (%s/%s)" % (vcf_record.REF,
                              vcf_record.samples[0]['AD'][0],
                              vcf_record.samples[0]['DP'])
        if len(vcf_record.ALT[0]) == 1:
            alt = "%s (%s/%s)" % (vcf_record.ALT[0],
                                  vcf_record.samples[0]['AD'][1],
                                  vcf_record.samples[0]['DP'])
        else:
            alt = "%sbp (%s/%s)" % (len(vcf_record.ALT[0]),
                                    vcf_record.samples[0]['AD'][1],
                                    vcf_record.samples[0]['DP'])
        filter_status = []

        # if any of the test failed, set PASS as failed
        if len(vcf_record.FILTER) != 0:
            vcf_record.FILTER.append('PASS')

        for filter_name in vcf_reader.filters:
            if filter_name in vcf_record.FILTER:
                if filter_name == 'PASS':
                    filter_status.append('NO')
                else:
                    filter_status.append('-')
            else:
                if filter_name == 'PASS':
                    filter_status.append('YES')
                else:
                    filter_status.append('+')
        row = [contig_name,
               len(contig),
               position,
               ref,
               alt,
               variant_feature["mut_location"],
               variant_feature["mut_type"],
               variant_feature["orf_name"],
               variant_feature["gene"],
               orf_before,
               orf_after]

        row += list(filter_status)

        #  if comparison to assembled genome, add data about self mapping
        #  (IF A VARIANT IS ALSO IDENTIFIED IN THAT MAPPING, PROBABLY A FALSE POSITIVE)
        if 'assembled' in vcf_file:
            GT, PASS = check_reference_mapping_GT(merged_vcf_records,
                                                  contig_name,
                                                  position,
                                                  reference)
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

vcf_reader = vcf.Reader(codecs.open(vcf_file, 'r', 'latin-1'))
n_snps_record = len([vcf_record for vcf_record in vcf_reader if vcf_record.samples[0]['GT'] not in ['.', '0']])
print("n_snps_record", n_snps_record)
if n_snps_record > 200:
    snp_table = "too much snps"
else:
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
with open(report_file, "w") as fh:
    publish_file(
        source=io.StringIO(report_str),
        destination=fh,
        writer_name="html",
        settings_overrides={"stylesheet_path": ""},
    )
