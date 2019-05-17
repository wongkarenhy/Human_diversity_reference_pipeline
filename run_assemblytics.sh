#!/bin/bash
cores="$1"
workdir="$2"
assemblytics="$3"
mode="$4"

# Exit immediately upon error
set -e

cd "$workdir"/assemblytics

while read -r SAMPLE SEX FASTQ_DIR LONGRANGER_DIR ASSM_DIR BN_DIR ENZYME SUPERNOVA_VER ALT_NAME NUCMER_DIR POPULATION; do 
    
    if [ "$mode" = "PB" ]; then
        
        (bash "${assemblytics}/Assemblytics" "$NUCMER_DIR" \
                "$SAMPLE" 10000; 

        python2.7 "${assemblytics}/compute_overlapping_coords.py" \
            -d "$NUCMER_DIR" \
            -i "$SAMPLE".variants_between_alignments.bed \
            -o "$SAMPLE".variants_between_alignments_new.bed) &
    
    else
    
        for haplo in 2.1 2.2; do

            # Use Walfred's assemblytics script
            (bash "${assemblytics}/Assemblytics" "$NUCMER_DIR""$haplo".delta \
                "$SAMPLE"_"$haplo" 10000; 

            python2.7 "${assemblytics}/compute_overlapping_coords.py" \
                -d "$NUCMER_DIR""$haplo".delta \
                -i "$SAMPLE"_"$haplo".variants_between_alignments.bed \
                -o "$SAMPLE"_"$haplo".variants_between_alignments_new.bed) &

        done
    
    fi
    
    background=( $(jobs -p) )
    if (( ${#background[@]} == cores )); then
        wait
    fi

done < "${workdir}/${mode}_sample_metadata.txt"    
