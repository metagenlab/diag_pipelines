tail -n +2 ${pipeline_folder}/data/validation_datasets/datasets.tsv | while read bioproject genus species therest
do
   mkdir -p ${bioproject}_${genus}-${species}/links/
   cd ${bioproject}_${genus}-${species}
   esearch -db sra -query "${bioproject}[BIOPROJECT] AND \"${genus} ${species}\"[ORGANISM]" < /dev/null | efetch -db sra -format runinfo | sed "s/,/\t/g" > ${bioproject}_${genus}-${species}.tsv
   ln -s ${main}/core_genomes/${genus}_${species}/ core_genome
   snakemake --snakefile ${pipeline_folder}/workflows/validation.rules --config sra_samples=${bioproject}_${genus}-${species}.tsv species="${genus} ${species}" absolute_path_of_pipeline_folder=${pipeline_folder} logging_folder="logging/" --configfile ${pipeline_folder}/data/validation_datasets/config.yaml -n
   cd ..
done 
