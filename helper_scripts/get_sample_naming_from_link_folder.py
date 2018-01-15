import os
import re


def get_read_naming_patterns_general(directory):
    result = []
    extension= {}
    for fname in os.listdir(directory):
        if fname.endswith("fastq.gz") or fname.endswith("fq.gz") or fname.endswith("fastq") or fname.endswith("fq"):
            print(fname)
            ext = re.search(r'_(R)?(1|2)(\.|_).*', fname)
            print(ext)
            if ext is None:
                ext = re.search(r'f(?:ast)?q(?:\.gz)?', fname)
                samp = re.sub("\.$", "", re.search(r'^([^\.]*)\.*', fname).group(0))
                if samp in extension.keys():
                    if ext.group(0).endswith(".gz"):
                        extension[samp]=[ext.group(0)]
                else:
                    extension[samp]=[ext.group(0)]
            else:
                read = re.compile(re.search(r'_(R)?(1|2)(\.|_)', fname).group(0).replace("_", "").replace(".", ""))
                samp = re.sub(r'_(R)?(1|2)$', '', re.search(r'.*_(R)?(1|2)', fname).group(0))
                extension.setdefault(samp, [])
                match = list(filter(read.match, extension[samp]))
                if len(match) == 1:
                    if match[0].endswith("fastq") or match[0].endswith("fq"):
                        extension[samp].remove(match[0])
                        extension[samp].append(re.sub("^_", "", ext.group(0)))
                else:
                    extension[samp].append(re.sub("^_", "", ext.group(0)))
    return(extension)


def get_read_naming_patterns_molis(directory):
    result = []
    extension= {}
    for fname in os.listdir(directory):
        if fname.endswith("fastq.gz") or fname.endswith("fq.gz") or fname.endswith("fastq") or fname.endswith("fq"):
            ext = re.search(r'_(R)?(1|2)(\.|_).*', fname)
            if ext is None:
                ext = re.search(r'f(?:ast)?q(?:\.gz)?', fname)
                samp = re.sub("\.$", "", re.search(r'^([^\.]*)\.*', fname).group(0))
                if samp in extension.keys():
                    if ext.group(0).endswith(".gz"):
                        extension[samp]=[ext.group(0)]
                else:
                    extension[samp]=[ext.group(0)]
            else:
                read = re.compile(re.search(r'_R?(1|2)(\.|_)', fname).group(0).replace("_", "").replace(".", ""))
                samp = re.sub(r'_R?(1|2)$', '', re.search(r'.*_R?(1|2)', fname).group(0)).split("_")[0]
                extension.setdefault(samp, [])
                match = list(filter(read.match, extension[samp]))
                if len(match) == 1:
                    if match[0].endswith("fastq") or match[0].endswith("fq"):
                        extension[samp].remove(match[0])
                        extension[samp].append(re.sub("^_", "", ext.group(0)))
                else:
                    extension[samp].append(re.sub("^_", "", ext.group(0)))
            
    return(extension)
