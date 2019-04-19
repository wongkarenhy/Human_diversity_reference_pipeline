#!/bin/bash
WORKDIR="$1"
METADATA="$2"

if [ -f "$WORKDIR"/discovery/supernova_idx.txt ]; then
    rm "$WORKDIR"/discovery/supernova_idx.txt
fi

while read -r SAMPLE SEX FASTQ_DIR LONGRANGER_DIR ASSM_DIR BN_DIR ENZYME SUPERNOVA_VER ALT_NAME NUCMER_DIR POPULATION; do
    
    path_sample_name=$(sed 's:.*/::' <<< "$ASSM_DIR" | cut -d_ -f1)
    
    for haplo in 2.1 2.2; do
    
        awk -v s="$SAMPLE" -v h="$haplo" 'BEGIN{FS=" "}{print $1,$2,s,h}' \
            "$ASSM_DIR"/"$path_sample_name"_pseudohap"$haplo".fasta.fai >> "$WORKDIR"/discovery/supernova_idx.txt
            
    done
    
done < "$METADATA"
