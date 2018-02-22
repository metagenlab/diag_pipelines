FROM cashalow/snakemake

WORKDIR ${pipeline_folder}

RUN git pull 

RUN mkdir -p ${main}/validation/Listeria_monocytogenes/

WORKDIR ${main}/validation/Listeria_monocytogenes/

RUN ln -s /home/pipeline_user/core_genomes/Listeria_monocytogenes/ core_genome

RUN cp /home/pipeline_user/snakemake_pipeline/data/validation_datasets/L_monocytogenes_ASM_NGS_2015/* .

RUN mkdir links
       
RUN snakemake --snakefile ${pipeline_folder}/workflows/typing.rules --use-conda --conda-prefix ${conda_folder} --configfile config.yaml --config ref_ids_for_mapping="" -j 4 typing/freebayes_joint_genotyping/core_ridom/264498/bwa/distances_snp.tsv 

RUN snakemake --snakefile targets.rules --use-conda --conda-prefix ${conda_folder} --configfile config.yaml --config ref_ids_for_mapping="" -j 4 -R call_variants_freebayes_vcf_all_samples

RUN chown pipeline_user -R ${main}

USER pipeline_user
