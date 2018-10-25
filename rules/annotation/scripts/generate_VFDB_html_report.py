
import pandas

blast_result = snakemake.input[0]
vfdb_annotation = snakemake.input[1]
sample = snakemake.params[0]
report_file = snakemake.output[0]

# ORF_ID
# 0 matching_sequence
# 1 percentage_identity
# 2 alignment_length
# 3 mismatch_number
# 4 gap_number
# 5 alignment_start_on_virulence_factor
# 6 alignment_end_on_virulence_factor
# 7 alignment_start_on_matching_sequence
# 8 aligment_end_on_matching_sequence
# 9 e-value
# 10 bitscore
# 11 amino_acid_sequence_of_matching_sequence
# 12 query_coverage
# 13 algorithm

# seq_id
# 0 vf_id
# 1 gene
# 2 bacteria
# 3 vf_description
# 4 product
# 5 gbk_accession


STYLE = """
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css" integrity="sha384-BVYiiSIFeK1dGmJRAkycuHAHRg32OmUcww7on3RYdg4Va+PmSTsz/K68vbdEjh4u" crossorigin="anonymous">
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap-theme.min.css" integrity="sha384-rHyoN1iRsVXV4nD0JutlnGaslCJuC7uwjduW9SVrLvRYooPp2bWYgmgJQIXwl/Sp" crossorigin="anonymous">
    <link rel="stylesheet" href="https://cdn.datatables.net/1.10.16/css/jquery.dataTables.min.css">     
    """

SCRIPT = """ 
    <script src="https://code.jquery.com/jquery-1.11.0.min.js"></script>      
    <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js" integrity="sha384-Tc5IQib027qvyjSMfHjOMaLkfuWVxZxUPnCJA7l2mCWNIpG9mGCD8wGNIcPD7Txa" crossorigin="anonymous"></script>
    <script type="text/javascript" src="https://cdn.datatables.net/1.10.16/js/jquery.dataTables.min.js"></script>

    <script>
    $(document).ready(function() {
        $('#VF_table').DataTable( {
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

HTML = '''
<!DOCTYPE html>
<html>
    <head>  
    %s
    %s
    </head>  

    <body>
    <div id='vftable' style="max-width:90%%; padding-left:35px;">
        <h1>Virulence Factors of sample: %s</h1>
        %s
    </div>    
    </body>
<html>
'''

def parse_VF_annotation(vfdb_annotation):
    vf_accession2annotation = {}
    with open(vfdb_annotation, 'r') as f:
        for row in f:
            data = row.split('\t')
            vf_accession2annotation[data[0]] = data[1:]
    return vf_accession2annotation

def prepare_dataframe(blast_result, vf_accession2annotation):
    vf_gene_template = '<a href="http://www.mgc.ac.cn/cgi-bin/VFs/gene.cgi?GeneID=%s">%s</a>'
    vf_template = '<a href=http://www.mgc.ac.cn/cgi-bin/VFs/vfs.cgi?VFID=%s#%s>%s</a> (%s)'
    VFDB_link = '<a href=https://www.vfdb.chlamdb.ch/VF_VFDB/%s>Link</a>'

    blast_result_table = pandas.read_csv(blast_result, delimiter="\t", header=0, index_col=0)

    vf_data = []
    header = ['ORF ID', "Hit", "Virulence Factor", "Link Conservation", "Product", "Organism", "Identity (%)", "e-value", "Bitscore"]
    for protein_id, blast_row in blast_result_table.iterrows():
        vf_accession  = blast_row["matching_sequence"]
        percentage_identity =  blast_row["percentage_identity"]
        bitscore =  blast_row["bitscore"]
        evalue = blast_row["e-value"]

        annotation = vf_accession2annotation[vf_accession]
        vf_id = annotation[0]
        vf_description = annotation[3]
        product = annotation[4]
        bacteria =annotation[2]

        vf_data.append([protein_id,
                        vf_gene_template % (vf_accession, vf_accession),
                        vf_template % (vf_id, vf_id,vf_id, vf_description),
                        VFDB_link % vf_accession,
                        product,
                        bacteria,
                        percentage_identity,
                        evalue,
                        bitscore])

    df = pandas.DataFrame(vf_data, columns=header)
    # cell content is truncated if colwidth not set to -1
    pandas.set_option('display.max_colwidth', -1)

    df_str = df.to_html(
        index=False,
        bold_rows=False,
        classes=["display", "dataTable"],
        table_id="VF_table",
        escape=False,
        border=0)
    return df_str

vf_accession2annotation = parse_VF_annotation(vfdb_annotation)

with open(report_file, 'w') as f:
    f.write(HTML % (SCRIPT,
                    STYLE,
                    sample,
                    prepare_dataframe(blast_result,
                                      vf_accession2annotation)))