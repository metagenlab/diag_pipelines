
# from snakemake.utils import report

# with open(input[0]) as vcf:
#    n_calls = sum(1 for l in vcf if not l.startswith("#"))
# n_samples = list(read_naming.keys()
import pandas

qualimap_report = snakemake.input["qualimap_report"]
ete_figure = snakemake.input["ete_figure"]
ete_figure_counts = snakemake.input["ete_figure_counts"]
virulence_reports = snakemake.input["virulence_reports"]
ordered_samples = snakemake.params["samples"]

output_file = snakemake.output[0]
blast_files = [pandas.read_csv(name, delimiter='\t') for name in snakemake.input["blast_results"]]

sample2n_VFs = {}
for n, sample in enumerate(ordered_samples):
    sample2n_VFs[sample] = len(blast_files[n])

print(sample2n_VFs)

sample2mlst={}


report_template = '''

<!DOCTYPE html>


<html>
    <head>        

    <script src="https://code.jquery.com/jquery-1.11.0.min.js"></script>

    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css" integrity="sha384-BVYiiSIFeK1dGmJRAkycuHAHRg32OmUcww7on3RYdg4Va+PmSTsz/K68vbdEjh4u" crossorigin="anonymous">
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap-theme.min.css" integrity="sha384-rHyoN1iRsVXV4nD0JutlnGaslCJuC7uwjduW9SVrLvRYooPp2bWYgmgJQIXwl/Sp" crossorigin="anonymous">

    <!-- Latest compiled and minified JavaScript -->
    <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js" integrity="sha384-Tc5IQib027qvyjSMfHjOMaLkfuWVxZxUPnCJA7l2mCWNIpG9mGCD8wGNIcPD7Txa" crossorigin="anonymous"></script>

    <link rel="stylesheet" href="https://cdn.datatables.net/1.10.16/css/jquery.dataTables.min.css">
    <script type="text/javascript" src="https://cdn.datatables.net/1.10.16/js/jquery.dataTables.min.js"></script

    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/c3/0.6.7/c3.min.css">
    <script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/c3/0.6.7/c3.min.js"></script>

    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/c3/0.6.7/c3.min.css">
    <script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/d3/5.7.0/d3.min.js"></script>
    
    <script type="text/javascript" src="https://cdnjs.cloudflare.com/ajax/libs/Chart.js/1.0.2/Chart.min.js"></script>

    <style>

/* make sidebar nav vertical */
@media (min-width: 768px){
  .affix-content .container {
    width: 700px;
  }

  html,body{
    background-color: #f8f8f8;
    height: 100%%;
    overflow: hidden;
  }
    .affix-content .container .page-header{
    margin-top: 0;
  }
  .sidebar-nav{
        position:fixed;
        width:100%%;
  }
  .affix-sidebar{
    padding-right:0;
    font-size:small;
    padding-left: 0;
  }
  .affix-row, .affix-container, .affix-content{
    height: 100%%;
    margin-left: 0;
    margin-right: 0;
  }
  .affix-content{
    background-color:white;
  }
  .sidebar-nav .navbar .navbar-collapse {
    padding: 0;
    max-height: none;
  }
  .sidebar-nav .navbar{
    border-radius:0;
    margin-bottom:0;
    border:0;
  }
  .sidebar-nav .navbar ul {
    float: none;
    display: block;
  }
  .sidebar-nav .navbar li {
    float: none;
    display: block;
  }
  .sidebar-nav .navbar li a {
    padding-top: 12px;
    padding-bottom: 12px;
  }
}

@media (min-width: 769px){
  .affix-content .container {
    width: 600px;
  }
    .affix-content .container .page-header{
    margin-top: 0;
  }
}

@media (min-width: 992px){
  .affix-content .container {
  width: 900px;
  }
    .affix-content .container .page-header{
    margin-top: 0;
  }
}

@media (min-width: 1220px){
  .affix-row{
    overflow: hidden;
  }

  .affix-content{
    overflow: auto;
  }

  .affix-content .container {
    width: 1000px;
  }

  .affix-content .container .page-header{
    margin-top: 0;
  }
  .affix-content{
    padding-right: 30px;
    padding-left: 30px;
  }
  .affix-title{
    border-bottom: 1px solid #ecf0f1;
    padding-bottom:10px;
  }
  .navbar-nav {
    margin: 0;
  }
  .navbar-collapse{
    padding: 0;
  }
  .sidebar-nav .navbar li a:hover {
    background-color: #428bca;
    color: white;
  }
  .sidebar-nav .navbar li a > .caret {
    margin-top: 8px;
  }
}


* {box-sizing: border-box;}

body {
  margin: 0;
  font-family: Arial, Helvetica, sans-serif;
}

.topnav {
  overflow: hidden;
  background-color: #e9e9e9;
  position: top;
}

.topnav a {

  float: left;
  display: block;
  color: black;
  text-align: center;
  padding: 14px 16px;
  text-decoration: none;
  font-size: 17px;
}

.pageheader {
  background-color: #e9e9e9;
  color: red;
  padding: 17px 16px;
  font-size: 22px;
}


.topnav a:hover {
  background-color: #ddd;
  color: black;
}

.topnav a.active {
  background-color: #2196F3;
  color: white;
}

.topnav .search-container {
  float: right;
}

.topnav input[type=text] {
  padding: 6px;
  margin-top: 8px;
  font-size: 17px;
  border: none;
}

.topnav .search-container button {
  float: right;
  padding: 6px 10px;
  margin-top: 8px;
  margin-right: 16px;
  background: #ddd;
  font-size: 17px;
  border: none;
  cursor: pointer;
}

.topnav .search-container button:hover {
  background: #ccc;
}

@media screen and (max-width: 600px) {
  .topnav .search-container {
    float: none;
  }
  .topnav a, .topnav input[type=text], .topnav .search-container button {
    float: none;
    display: block;
    text-align: left;
    width: 100%%;
    margin: 0;
    padding: 14px;
  }
  .topnav input[type=text] {
    border: 1px solid #ccc;
  }
}
body {
    padding-bottom: 50px;
}
    </style>
</head>

<div class="topnav">
  
</div>


<body style="position: absolute;">

    <div class="col-sm-2 col-md-2 affix-sidebar" style="padding-top: 0px;">
      <div class="sidebar-nav" >
        <div class="navbar navbar-default" role="navigation">
            <div class="navbar-collapse collapse sidebar-navbar-collapse">
              <ul class="nav navbar-nav" id="sidenav01">
                    <li class="active">

                    </li>
                    <li class="pageheader"><a>Virulence/Typing <br> 
                    Report</a></li>
                    <li><a href="#quality">1. QC</a></li>
                    <li><a href="#phylogeny">2. MLST</a></li>
                    <li><a href="-">3. Antibiotic Resistances</a></li>
                    <li><a href="-">4. Virulence</a></li>
                    <li><a href="-">5. Phylogeny</a></li>
                    <li><a href="-">6. Minimum spanning tree</a></li>
              </ul>
            </div><!--/.nav-collapse -->
         </div>
      </div>
    </div>

    <div class="col-sm-10 col-md-10 affix-content">
            <h1 id="quality">1. Quality control</h1>
            <ul>
                <li><a href="%s">MULTIQC</a></li>
            </ul>
            <h1 id="phylogeny">2. MLST</h1>
            <div>
                <img style="width:90%%" src="%s" align="top">
            </div>
            <h1 id="phylogeny">3. Antibiotic Resistances</h1>
            <h1 id="phylogeny">4. Virulence factors (VFDB)</h1>
            
            <h3>4.1 Overview</h3>
            <img style="width:50%%" src="%s" align="top">
            <h3>4.2 Details</h3>
            %s
            <h1 id="phylogeny">5. Phylogeny</h1>
            
            <img style="width:90%%" src="%s" align="top">
            
            <h1 id="phylogeny">6. Minimum spanning tree</h1>
    </div>




</body>
<script>


var chart = c3.generate({
    size: {
        height: 400,
        width: 800
    },
    legend: {
        show: false
    },
    grid: {
        x: {
            show: false
        },
        y: {
            show: true
        }
    },
    data: {
        columns: [
            [ "# of VFs",3,1,2,3,1,1,1,2,1,2,1,1,1,3,2,1,1,2,1,1,2,2,1,1,2,2,2,1,1,3,1,1,3,3,6,4,3,10,20,28,75]
        ],
        type: 'bar',
        onclick: function (e) {
            console.log(e);
            console.log(this.internal.config.axis_x_categories[e.index]);
        },
    },
    bar: {
        width: {
            ratio: 0.5 // this makes bar width 50%% of length between ticks
        }
        // or
        //width: 100 // this makes bar width 100px
    },
    axis: {
      rotated: false,
      x: {
        type: 'categorized',
        categories: ["1","4","5","9","13","15","17","18","31","36","54","58","59","60","61","62","63","76","107","121","123","133","136","150","151","152","158","160","163","166","168","169","170","172","174","175","176","177","178","179","180"],
        tick: {
            rotate: 75,
            multiline: false
        },
        label: {
            text: 'Number of genomes',
            position: 'outer-center'
            // inner-right : default
            // inner-center
            // inner-left
            // outer-right
            // outer-center
            // outer-left
        }
     },
      y: {
        label: {
            text: 'VF count',
            position: 'outer-middle'
            // inner-right : default
            // inner-center
            // inner-left
            // outer-right
            // outer-center
            // outer-left
        }
     }

    }
});

</script>

        '''

virulence_section = '''
                            <table class="display dataTable" id="VF_table">
                                <thead>
                                    <tr>
                                      <th scope="col">Strain id</th></th>
                                      <th scope="col">Number of VFs</th></th>
                                      <th scope="col">Virulence Report</th>
                                </thead>
                                <tbody>
                                    %s
                                </tbody>
                            </table>        
                            '''

virulence_row = '''
                        <tr>
                            <td>%s</td>
                            <td>%s</td>
                            <td><a href="%s">%s</a></td>

                        </tr>        
                        '''
rows = ''
for report in virulence_reports:
    print(report.split('/'))
    rows += virulence_row % (report.split('/')[1], sample2n_VFs[report.split('/')[1]],report, report)

with open(output_file, 'w') as f:
    f.write(report_template % (qualimap_report,
                               ete_figure,
                               ete_figure_counts,
                               virulence_section % rows,
                               ete_figure,))
f.close()
