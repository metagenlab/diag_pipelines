
import pandas
import re

rgi_tsv_output = snakemake.input[0]
rgi_ontology = snakemake.input[1]
gene_depth_file = snakemake.input[2]
contig_gc_depth_file = snakemake.input[3]
samtools_depth = snakemake.input[4]
sample = snakemake.params[0]

report_file = snakemake.output[0]

# parse rgi output file
orf_id2data = {}
contig2resistances = {}
with open(rgi_tsv_output, 'r') as f:
    for row in f:
        data = row.split('\t')
        contig = '_'.join(data[0].split('_')[0:-1])
        if contig not in contig2resistances:
            contig2resistances[contig] = [re.sub("'", "", data[8])]
        else:
            contig2resistances[contig].append(re.sub("'", "", data[8]))
        orf_id2data[data[0]] = data[1:len(data)]
print("contig2resistances")
print(contig2resistances)

# calculate rgi hit(s) sequencing depth based on position of the CDS in contigs
def parse_smatools_depth(samtools_depth):
    import pandas

    with open(samtools_depth, 'r') as f:
        table = pandas.read_csv(f, sep='\t', header=None, index_col=0)
    return table

samtools_dataframe = parse_smatools_depth(samtools_depth)

# get contig depth and GC
contig2gc_content = pandas.read_csv(contig_gc_depth_file,
                                      delimiter='\t',
                                      header=0).set_index("contig").to_dict()["gc_content"]

contig2median_depth = pandas.read_csv(contig_gc_depth_file,
                                      delimiter='\t',
                                      header=0).set_index("contig").to_dict()["median_depth"]

contig2size = pandas.read_csv(contig_gc_depth_file,
                                      delimiter='\t',
                                      header=0).set_index("contig").to_dict()["contig_size"]
# prepare bubble plot GC vs Coverage

dataset_template = '''
        {
        label: '%s',
        data: [
            %s
        ],
        backgroundColor:"%s",
        hoverBackgroundColor: "%s"
        }
'''

data_template = '''{
            x: %s,
            y: %s,
            r: Math.log(%s),
            l: "%s",
          }'''

bubble_template = '''
var data = {
    datasets: [
      %s,
      %s
      ]
  };

var options = {
         legend: {
            display: false
         },
          type: 'bubble',
          data: data,
          options: {
            tooltips: {
              enabled: true,
              callbacks: {
                label: function(tooltipItem, data) {
                  var label = data.datasets[tooltipItem.datasetIndex].data[tooltipItem.index].l;
                  return label;
                }
         
            }
            },
          title:{
              display: true,
              text:'Depth vs GC plot'
          },
            scales: {
                xAxes: [{
                        display: true,
                        scaleLabel: {
                            display: true,
                            labelString: 'GC(%%)'
                        },
                                                    ticks: {
                            beginAtZero: true,
                            steps: 10,
                            stepValue: 10,
                            max: 100
                        }
                    }],
                yAxes: [{
                        display: true,
                        scaleLabel: {
                            display: true,
                            labelString: 'Depth'
                        },
                                                    ticks: {
                            beginAtZero: true,
                        }
                    }]
            }, 
        // Container for pan options
        pan: {
            // Boolean to enable panning
            enabled: true,

            // Panning directions. Remove the appropriate direction to disable 
            // Eg. 'y' would only allow panning in the y direction
            mode: 'xy'
        },

        // Container for zoom options
        zoom: {
            // Boolean to enable zooming
            enabled: true,

            // Zooming directions. Remove the appropriate direction to disable 
            // Eg. 'y' would only allow zooming in the y direction
            mode: 'x',
        }

              
        }
}
var ctx = document.getElementById('chartJSContainer').getContext('2d');
new Chart(ctx, options);
'''

no_resistance_list = []
resistance_list = []



for contig in contig2gc_content:
    depth = contig2median_depth[contig]
    gc = contig2gc_content[contig]
    if contig in contig2resistances:
        resistance_list.append(data_template % (gc,
                                                depth,
                                                contig2size[contig], # contig size
                                                "%s (%sbp): %s" % (contig,
                                                                   contig2size[contig],
                                                                   ','.join(contig2resistances[contig]))))
    else:
        no_resistance_list.append(data_template % (gc,
                                                depth,
                                                contig2size[contig], # contig size
                                                "%s (%sbp)" % (contig,
                                                               contig2size[contig])))

buuble_chart_code = bubble_template % (dataset_template % ('No resistances',
                                                           ','.join(no_resistance_list),
                                                           "rgba(22, 99, 132, 0.2)",
                                                           "rgba(22, 99, 132, 0.2)"),
                                       dataset_template % ('Resistances',
                                                           ','.join(resistance_list),
                                                           '#ff6384',
                                                           '#ff6384'
                                                           ))

# parse rgi ontology



template = '''
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

    <script type="text/javascript" src="https://github.com/chartjs/Chart.js/releases/download/v2.7.1/Chart.bundle.min.js"></script>
<script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/chartjs-plugin-zoom/0.6.5/chartjs-plugin-zoom.js"></script>

    </head>
<body style="max-width:95%%; padding-left:55px;">
    <h1>Identification of resistance markers</h1>
    
    <h2>GC-coverage plot</h2>
        <div style="width:400px;">
        <canvas id="chartJSContainer" width="400" height="400"></canvas>       
        </div>
    
    <h2>Detailed table<h2>
        <div id='vftable' style="max-width:90%%; padding-left:55px;">
            <table class="display dataTable" id="res_table">
                <thead>
                    <tr>
                      <th scope="col">ARO</th></th>
                      <th scope="col">Name</th>
                      <th scope="col">Gene</th>
                      <th scope="col">Resistance Type</th>
                      <th scope="col">Variant</th>
                      <th scope="col">Organism</th>
                      <th scope="col">Identity (%%)</th>
                      <th scope="col">e-value</th>
                      <th scope="col">bit-score</th>
                </thead>
                <tbody>
    
                </tbody>
            </table>
        </div>
    <h2>Circos<h2>
    

</body>
<script type="text/javascript">
%s

$(document).ready(function() {
    $('#res_table').DataTable( {
        dom: 'Bfrtip',
        "pageLength": 20,
        "searching": true,
        "bLengthChange": false,
        "paging":   true,
        "info": false
    } );
} );


</script>

</html>

'''

row_template = '''
    <tr>
        <td>%s</td>
        <td>%s</td>
        <td>%s</td>
        <td>%s</td>
        <td>%s</td>
        <td>%s</td>
        <td>%s</td>
        <td>%s</td>
    </tr>
'''

with open(report_file, 'w') as f:
    f.write(template % (buuble_chart_code))