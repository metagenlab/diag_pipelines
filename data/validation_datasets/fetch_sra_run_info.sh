tail -n +2 ${pipeline_folder}/data/validation_datasets/datasets.tsv | while read bioproject genus species therest
do
    mkdir -p ${bioproject}_${genus}-${species}/links/
    echo ${main}/core_genomes/${genus}_${species%%_*}
    cd ${bioproject}_${genus}-${species}
    ln -s -f ${main}/core_genomes/${genus}_${species%%_*} core_genome
    esearch -db sra -query "${bioproject}[BIOPROJECT] AND \"${genus} ${species}\"[ORGANISM]" < /dev/null | efetch -db sra -format runinfo | sed "s/,/\t/g" | head -n 1 > ${bioproject}_${genus}-${species}.tsv
    esearch -db sra -query "${bioproject}[BIOPROJECT] AND \"${genus} ${species}\"[ORGANISM]" < /dev/null | efetch -db sra -format runinfo |  grep -v "Run,ReleaseDate," | sed "s/,/\t/g" | grep -v "^$" >> ${bioproject}_${genus}-${species}.tsv
    snakemake --snakefile ${pipeline_folder}/workflows/typing.rules --config sra_samples=${bioproject}_${genus}-${species}.tsv species="${genus} ${species}" --configfile ${pipeline_folder}/data/validation_datasets/config.yaml -n typing/freebayes_joint_genotyping/cgMLST/bwa/distances_in_snp_mst_no_st.svg typing/gatk_gvcfs/cgMLST/bwa/distances_in_snp_mst_no_st.svg phylogeny/freebayes_joint_genotyping/cgMLST/bwa/phylogeny_no_st.svg

    cd ..
done 
