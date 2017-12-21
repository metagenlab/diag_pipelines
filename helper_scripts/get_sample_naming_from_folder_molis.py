import os
import re

def get_read_naming_patterns(directory):
    result = []
    extension= {}
    for fname in os.listdir(directory):
        if fname.endswith("fastq.gz") or fname.endswith("fq.gz") or fname.endswith("fastq") or fname.endswith("fq"):
            ext = re.search(r'_.*_R(1|2)(\.|_).*', fname)
            read = re.compile(re.search(r'_R(1|2)(\.|_)', fname).group(0).replace("_", "").replace(".", ""))
            samp = re.sub(r'_R(1|2)$', '', re.search(r'.*_R(1|2)', fname).group(0)).split("_")[0]
            extension.setdefault(samp, [])
            match = list(filter(read.match, extension[samp]))
            if len(match) == 1:
                if match[0].endswith("fastq") or match[0].endswith("fq"):
                    extension[samp].remove(match[0])
                    extension[samp].append(re.sub("^_", "", ext.group(0)))
            else:
                extension[samp].append(re.sub("^_", "", ext.group(0)))
    return(extension)



