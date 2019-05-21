#!/bin/bash
WORKDIR="$1"

if [ -f "$WORKDIR"/discovery/tmp_idx.txt ]; then
    rm "$WORKDIR"/discovery/tmp_idx.txt
fi

while read -r SAMPLE SEX FASTQ_DIR LONGRANGER_DIR ASSM_DIR BN_DIR ENZYME SUPERNOVA_VER ALT_NAME NUCMER_DIR POPULATION SRC; do
    
    if [ "$SRC" = "10X" ]; then
    
        path_sample_name=$(sed 's:.*/::' <<< "$ASSM_DIR" | cut -d_ -f1)

        for haplo in 2.1 2.2; do

            if [ ! -f "$ASSM_DIR"/"$path_sample_name"_pseudohap"$haplo".fasta.fai ]; then
                samtools faidx "$ASSM_DIR"/"$path_sample_name"_pseudohap"$haplo".fasta
            fi
        
            awk -v s="$SAMPLE" -v h="$haplo" 'BEGIN{FS=" "}{print $1,$2,s,h}' \
                "$ASSM_DIR"/"$path_sample_name"_pseudohap"$haplo".fasta.fai >> "$WORKDIR"/discovery/tmp_idx.txt

        done
    
    elif [ "$SRC" = "PB" ]; then
        
        if [ ! -f "${ASSM_DIR}.fai" ]; then
            samtools faidx ${ASSM_DIR}
        fi
        
        awk -v s="$SAMPLE" 'BEGIN{FS=" "}{print $1,$2,s,"unphased"}' \
            "${ASSM_DIR}.fai" >> "$WORKDIR"/discovery/tmp_idx.txt
    
    else
    
        echo "Error: unknwon source!"
        exit 1
        
    fi
    
done < "${WORKDIR}/TMP_sample_metadata.txt"
