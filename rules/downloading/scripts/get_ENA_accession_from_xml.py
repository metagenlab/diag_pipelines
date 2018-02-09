from lxml import etree
import glob
import re


for i in glob.glob(str(snakemake.input).replace("/xmls_log.txt", "/xmls/")+"/*"):
    id = re.sub(".*/", "", i).replace(".xml", "")
    tree = etree.parse(i)
    xpath_dbref_EMBL = "/u:uniprot/u:entry/u:dbReference[@type='EMBL']"
    xpath_seq_ID = "u:property[@type='protein sequence ID' and ../u:property[@type='molecule type' and @value='Genomic_DNA']]"
    ena_id=""
    print(id)
    tree_req_EMBL=tree.xpath(xpath_dbref_EMBL, namespaces = {'u':'http://uniprot.org/uniprot'})
    if len(tree_req_EMBL):
        print(len(tree_req_EMBL))
        for item in tree_req_EMBL:
            if len(item.xpath(xpath_seq_ID, namespaces = {'u':'http://uniprot.org/uniprot'})):
                for item2 in item.xpath(xpath_seq_ID, namespaces = {'u':'http://uniprot.org/uniprot'}):
                    ena_id=item2.attrib["value"]
    else:
        raise ValueError("Missing EMBL entry in {0}".format(i))
    if not ena_id:
        raise ValueError("Missing ENA ID in {0}".format(i))
    else:
        with open(snakemake.output[0], "a") as f_ena:
            f_ena.write(id+"\t"+ena_id+"\n")



            




