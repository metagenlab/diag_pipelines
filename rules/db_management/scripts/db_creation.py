import mysql.connector

cnx = mysql.connector.connect(option_files=snakemake.params["conf"], option_groups=snakemake.params["db"])
cnx.get_warnings = True
cursor = cnx.cursor()

cmds={}
cmds['samples']="""
CREATE TABLE samples(Specimen varchar(255) NOT NULL, Patient varchar(255) NOT NULL, Bio_Project varchar(255)  NOT NULL, Material  varchar(255) NOT NULL, Specimen_Collected_Date Date NOT NULL, Bio_Sample varchar(255) NOT NULL, Lineage varchar(255) NOT NULL, Octal_Spoligotype varchar(255) NOT NULL, PRIMARY KEY (Specimen));
"""


cmds['sras']="""
CREATE TABLE sras(specimen varchar(255) NOT NULL, Sequence_Read_Archive varchar(255) NOT NULL, PRIMARY KEY(specimen, Sequence_Read_Archive));
"""

cmds["dst"]="""
CREATE TABLE dst(specimen varchar(255) NOT NULL, test varchar(255) NOT NULL, Date date, antibio varchar(255) NOT NULL, phenotype varchar(255) NOT NULL, PRIMARY KEY (specimen, test, Date, antibio));
"""

cmds["corres"]="""
CREATE TABLE corres(specimen varchar(255) NOT NULL, file varchar(255), PRIMARY KEY(specimen));
"""

cmds["resistance_associated_mutations"]="""
CREATE TABLE resistance_associated_mutations(specimen varchar(255) NOT NULL, software varchar(255) NOT NULL, gene varchar(255) NOT NULL, position int NOT NULL, ref_aa varchar(255) NOT NULL, mut_aa varchar(255) NOT NULL, PRIMARY KEY(specimen, software, gene, position));
"""

cmds["phenotype_prediction"]="""
CREATE TABLE phenotype_prediction(specimen varchar(255) NOT NULL, software varchar(255) NOT NULL, antibiotic varchar(255) NOT NULL, PRIMARY KEY(specimen, software, antibiotic));
"""

cmds["antibiotic_resistance_conferring_genes_annotation"]="""
CREATE TABLE antibiotic_resistance_conferring_genes_annotation (gene varchar(255) NOT NULL, annotation_source varchar(255) NOT NULL, antibiotic varchar(255), PRIMARY KEY(gene, annotation_source, antibiotic));
"""

cmds["antibiotic_resistance_conferring_mutations_annotation"]="""
CREATE TABLE antiobiotic_resistance_conferring_mutations_annotation (gene varchar(255) NOT NULL, position int NOT NULL, ref_aa varchar(255) NOT NULL, mut_aa varchar(255), annotation_source varchar(255) NOT NULL, antibiotic varchar(255), PRIMARY KEY(gene, position, ref_aa, mut_aa, annotation_source, antibiotic));
"""

cmds["resistance_associated_genes"]="""
CREATE TABLE resistance_associated_genes(specimen varchar(255) NOT NULL, software varchar(255) NOT NULL, gene varchar(255) NOT NULL, PRIMARY KEY(specimen, software, gene));
"""

with open(snakemake.output[0], "w") as myfile:
    myfile.write("")

for i in cmds.keys():
    try:
        cursor.execute(cmds[i])
    except mysql.connector.errors.Error as err :
        with open(snakemake.output[0], "a") as f:
            f.write("Something went wrong: {}\n".format(err))
        
cnx.commit()
cnx.close()


