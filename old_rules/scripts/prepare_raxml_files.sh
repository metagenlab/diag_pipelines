#!/bin/sh
target_folder=$(dirname {input})
type=prot
model=LG
index=0
start=1
for i in $(ls ${{target_folder}}/*_oneline.fa)
do
    length=$(tail -n -1 ${{i}} | tr -d '\\n' | wc -m);
    echo "${{model}}, ${{type}}${{index}} = ${{start}}-$(( start + length - 1 ))" >> {output[1]}
    start=$(( start + length ));
    index=$(( index + 1 ));
done
paste -d'\\0' ${{target_folder}}/*_oneline.fa | sed "s/>\([0-9]\+\)>.*/>\\1/"  |  sed "/^\s*$/d" | sed "s/U/C/g" > {output[0]} #replaces Selenocysteine (U), which is not accepted by raxml, by Cysteine (C)
