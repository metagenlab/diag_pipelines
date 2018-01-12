import hashlib

with open(snakemake.output[0], "w") as sha256_file:
    with open(snakemake.input[0], "r", encoding='utf-8') as xmfa_file:
        for line in xmfa_file.readlines():
            if line.startswith("#"):
                locus_tag = line.replace("#", "").strip()
            elif line.startswith(">"):
                allele = line.replace(">", "").strip()
            elif line.strip():
                sha256 = hashlib.sha256(line.strip().encode()).hexdigest()
                sha256_file.write(locus_tag+"_"+allele+"\t"+sha256+"\n")
                
            
        
            
            
