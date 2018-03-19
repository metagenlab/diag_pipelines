tail -n +2 ${pipeline_folder}/data/validation_datasets/datasets.tsv | while read bioproject genus species therest
do
    mkdir -p ${bioproject}_${genus}-${species}/links/
    echo ${main}/core_genomes/${genus}_${species%%_*}
    cd ${bioproject}_${genus}-${species}
    ln -s -f ${main}/core_genomes/${genus}_${species%%_*} core_genome
    esearch -db sra -query "${bioproject}[BIOPROJECT] AND \"${genus} ${species}\"[ORGANISM]" < /dev/null | efetch -db sra -format runinfo | sed "s/,/\t/g" | head -n 1 > ${bioproject}_${genus}-${species}.tsv
    esearch -db sra -query "${bioproject}[BIOPROJECT] AND \"${genus} ${species}\"[ORGANISM]" < /dev/null | efetch -db sra -format runinfo |  grep -v "Run,ReleaseDate," | sed "s/,/\t/g" | grep -v "^$" >> ${bioproject}_${genus}-${species}.tsv
    snakemake --snakefile ${pipeline_folder}/workflows/validation.rules --config sra_samples=${bioproject}_${genus}-${species}.tsv species="${genus} ${species}" absolute_path_of_pipeline_folder=${pipeline_folder} logging_folder="logging/" --configfile ${pipeline_folder}/data/validation_datasets/config.yaml -n -q
    cd ..
done 
