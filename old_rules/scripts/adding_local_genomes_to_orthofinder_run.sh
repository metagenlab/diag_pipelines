#!/bin/sh

#This is the script that adds to a previous orthofinder run the proteomes of the newly sequenced strains.
#Orthofinder first runs in a termporary folder and if the run is successful this folder replaces the main result folder.


input_folder=$(dirname {input[1]})
storing_folder=$(dirname ${{input_folder}})/added/
if [ ! "$(find ${{input_folder}} -name *.faa)" ] #checking whether proteome files can be added
then
    echo 'No additional fasta files included for new orthofinder run'
    cp $(dirname {log})/tmp/id.txt {output[0]}
    cp $(dirname {log})/tmp/log.txt {output[2]}
    sed "s/://" $(dirname {log})/tmp/WorkingDirectory/SpeciesIDs.txt | sed "s/\\.faa//" | sed "s/^#.*//" > {output[1]}
    #we copy again the files that have been deleted by the snakemake rule
else
    rm -rf $(dirname {log})/tmp/
    cp -R $(dirname {input[0]}) $(dirname {log})/tmp/ 
    orthofinder -oa -M msa -b $(dirname {log})/tmp/WorkingDirectory/ -f ${{input_folder}} > {log}
    grep -no "Statistics_Overall_[0-9]\\+.csv" {log} | sed "s/.*\\([0-9]\\+\\).csv/\\1/" > $(dirname {log})/tmp/id.txt #looking for the id number of the completed run (it is incremented each time a run completes)
    if [ -s $(dirname {log})/tmp/id.txt ]
    then
        rm -rf $(dirname {input[0]})
        mv {log} $(dirname {log})/tmp/log.txt 
        cp -R $(dirname {log})/tmp/ $(dirname {input[0]}) #if the run is successful we replace the old orthofinder results folder
        sed "s/://" $(dirname {output[0]})/WorkingDirectory/SpeciesIDs.txt | sed "s/\\.faa//" | sed "s/^#.*//" | sed "s/\[/{{/" | sed "s/\]/}}/" > {output[1]} #we store the species name in a readable way
        mkdir -p ${{storing_folder}}
        mv ${{input_folder}}/*.faa ${{storing_folder}} #we move the proteome files that have been successfully added to the orthofinder run to a new folder
    fi
fi
