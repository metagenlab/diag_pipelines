import os
import re


def get_read_naming_patterns(directory):
    result = []
    extension= {}
    for fname in os.listdir(directory):
        if fname.endswith("fastq.gz") or fname.endswith("fq.gz") or fname.endswith("fastq") or fname.endswith("fq"):
            regex_str = '(_L0+[1-9]+)?_(R)?(1|2)(\.|_)'
            regex = re.compile(regex_str)
            ext = re.search(regex, fname)
            if ext is None:
                ext = re.search(r'f(?:ast)?q(?:\.gz)?', fname)
                samp = re.sub("\.$", "", re.search(r'^([^\.]*)\.*', fname).group(0))
                if samp in extension.keys():
                    if ext.group(0).endswith(".gz"):
                        extension[samp]=[ext.group(0)]
                else:
                    extension[samp]=[ext.group(0)]
            else:
                regex_after = re.compile(regex_str+".*")
                regex_before = re.compile(".*"+regex_str)
                read = re.compile(re.search(regex_after, fname).group(0))
                samp = re.sub(regex, '', re.search(regex_before, fname).group(0))
                extension.setdefault(samp, [])
                extension[samp].append(re.sub("^_", "", read.pattern))
    return(extension)

