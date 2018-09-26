
blast_result = snakemake.input[0]
vfdb_annotation = snakemake.input[1]
sample = snakemake.params[0]

vf_accession2annotation = {}
with open(vfdb_annotation, 'r') as f:
    for row in f:
        data = row.split('\t')
        vf_accession2annotation[data[0]] = data[1:]
protein_id2blast_hit = {}
with open(blast_result, 'r') as f:
    for row in f:
        data = row.split('\t')
        if data[0] in ['ORF_ID', 'virulence_factor_ID']:
            continue
        protein_id2blast_hit[data[0]] = data[1:]

report_file = snakemake.output[0]

template = ''''
<!DOCTYPE html>


<html>
    <head>        
        
        <script src="https://code.jquery.com/jquery-1.11.0.min.js"></script>
        
        <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css" integrity="sha384-BVYiiSIFeK1dGmJRAkycuHAHRg32OmUcww7on3RYdg4Va+PmSTsz/K68vbdEjh4u" crossorigin="anonymous">
        <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap-theme.min.css" integrity="sha384-rHyoN1iRsVXV4nD0JutlnGaslCJuC7uwjduW9SVrLvRYooPp2bWYgmgJQIXwl/Sp" crossorigin="anonymous">
        
        <!-- Latest compiled and minified JavaScript -->
        <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js" integrity="sha384-Tc5IQib027qvyjSMfHjOMaLkfuWVxZxUPnCJA7l2mCWNIpG9mGCD8wGNIcPD7Txa" crossorigin="anonymous"></script>
        
        <link rel="stylesheet" href="https://cdn.datatables.net/1.10.16/css/jquery.dataTables.min.css">
        <script type="text/javascript" src="https://cdn.datatables.net/1.10.16/js/jquery.dataTables.min.js"></script>
    
    </head>
<body>
    <div id='vftable' style="max-width:90%%; padding-left:35;">
        <h1>Virulence Factors of sample: %s</h1>
        <table class="display dataTable" id="VF_table">
            <thead>
                <tr>
                    <th scope="col">ORF ID</th></th>
                  <th scope="col">Hit</th>
                  <th scope="col">Virulence Factor</th>
                  <th scope="col">Link conservation</th>
                  <th scope="col">Product</th>
                  <th scope="col">Organism</th>
                  <th scope="col">Identity (%%)</th>
                  <th scope="col">e-value</th>
                  <th scope="col">bit-score</th>
            </thead>
            <tbody>
                %s
            </tbody>
        </table>
    </div>

</body>

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


'''

row_template = '''
    <tr>
        <td>%s</td>
        <td><a href="http://www.mgc.ac.cn/cgi-bin/VFs/gene.cgi?GeneID=%s">%s</a></td>
        <td><a href=http://www.mgc.ac.cn/cgi-bin/VFs/vfs.cgi?VFID=%s#%s>%s</a> (%s)</td>
        <td><a href=https://www.vfdb.chlamdb.ch/VF_VFDB/%s>Link</a></td>
        <td>%s</td>
        <td>%s</td>
        <td>%s</td>
        <td>%s</td>
        <td>%s</td>
    </tr>
'''

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

table_rows = ''
for blast_hit in protein_id2blast_hit:
    vf_accession = protein_id2blast_hit[blast_hit][0]
    annotation = vf_accession2annotation[vf_accession]
    table_rows+=row_template % (blast_hit,
                                protein_id2blast_hit[blast_hit][0],
                                protein_id2blast_hit[blast_hit][0],
                                annotation[0],
                                annotation[0],
                                annotation[0],
                                annotation[3],
                                protein_id2blast_hit[blast_hit][0],
                                annotation[4],
                                annotation[2],
                                protein_id2blast_hit[blast_hit][1],
                                protein_id2blast_hit[blast_hit][9],
                                protein_id2blast_hit[blast_hit][10],)

with open(report_file, 'w') as f:
    f.write(template % (sample,
                        table_rows))