
import pandas

blast_result = snakemake.input[0]
sample = snakemake.params[0]
report_file = snakemake.output[0]


# virulence_factor_ID
# description	matching_sequence
# percentage_identity
# alignment_length
# mismatch_number
# gap_number
# alignment_start_on_virulence_factor
# alignment_end_on_virulence_factor
# alignment_start_on_matching_sequence
# aligment_end_on_matching_sequence
# e-value	bitscore
# amino_acid_sequence_of_matching_sequence
# query_coverage
# algorithm
# gene
# gene_uniprot


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

def prepare_dataframe(blast_result):
    uniprot_template = '<a href="https://www.uniprot.org/uniprot/%s">%s</a>'

    blast_result_table = pandas.read_excel(blast_result, header=0, index_col=0)

    vf_data = []
    header = ['uniprot_accession', 'VF ID', "description", "percentage_identity", "alignment_length", "e-value", "query_coverage", "algorithm", "gene_uniprot"]
    for protein_id, blast_row in blast_result_table.iterrows():
        description =  blast_row["description"]
        percentage_identity =  blast_row["percentage_identity"]
        alignment_length = blast_row["alignment_length"]
        evalue = blast_row["e-value"]
        query_coverage = blast_row["query_coverage"]
        algorithm = blast_row["algorithm"]
        gene_uniprot = blast_row["gene_uniprot"]
        uniprot_accession = protein_id.split("_")[1]

        vf_data.append([uniprot_template % (uniprot_accession, uniprot_accession),
                        protein_id,
                        description,
                        percentage_identity,
                        alignment_length,
                        evalue,
                        query_coverage,
                        algorithm,
                        gene_uniprot])

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


with open(report_file, 'w') as f:
    f.write(HTML % (SCRIPT,
                    STYLE,
                    sample,
                    prepare_dataframe(blast_result)))
