#!/bin/bash
cores="$1"
workdir="$2"
assemblytics="$3"

# Exit immediately upon error
set -e

cd "$workdir"/assemblytics

while read -r SAMPLE SEX FASTQ_DIR LONGRANGER_DIR ASSM_DIR BN_DIR ENZYME SUPERNOVA_VER ALT_NAME NUCMER_DIR POPULATION SRC; do 
    
    if [ "$SRC" = "PB" ]; then
        
        (bash "${assemblytics}/Assemblytics" "$NUCMER_DIR" \
                "$SAMPLE" 10000; 

        python2.7 "${assemblytics}/compute_overlapping_coords.py" \
            -d "$NUCMER_DIR" \
            -i "$SAMPLE".variants_between_alignments.bed \
            -o "$SAMPLE".variants_between_alignments_new.bed;
         
        sed -e 's/$/\t\./' "$SAMPLE".variants_within_alignments.bed | \
            cat - "$SAMPLE".variants_between_alignments_new.bed | \
            grep -E 'Insertion|Repeat_expanion|Tandem_expanion' | \
            awk '$5>=10' > "$SAMPLE".filtered_variants.bed;
         
        echo [`date +"%Y-%m-%d %H:%M:%S"`] "run_assemblytics.sh: ${SAMPLE}_unphased DONE") &
    
    else
    
        for haplo in 2.1 2.2; do

            # Use Walfred's assemblytics script
            (bash "${assemblytics}/Assemblytics" "$NUCMER_DIR""$haplo".delta \
                "$SAMPLE"_"$haplo" 10000; 

            python2.7 "${assemblytics}/compute_overlapping_coords.py" \
                -d "$NUCMER_DIR""$haplo".delta \
                -i "$SAMPLE"_"$haplo".variants_between_alignments.bed \
                -o "$SAMPLE"_"$haplo".variants_between_alignments_new.bed;
            
            sed -e 's/$/\t\./' "$SAMPLE"_"$haplo".variants_within_alignments.bed | \
                cat - "$SAMPLE"_"$haplo".variants_between_alignments_new.bed | \
                grep -E 'Insertion|Repeat_expanion|Tandem_expanion' | \
                awk '$5>=10' > "$SAMPLE"_"$haplo".filtered_variants.bed;

            echo [`date +"%Y-%m-%d %H:%M:%S"`] "run_assemblytics.sh: ${SAMPLE}_${haplo} DONE") &

        done
    
    fi
    
    background=( $(jobs -p) )
    if (( ${#background[@]} == cores )); then
        wait
    fi

done < "${workdir}/TMP_sample_metadata.txt"    
        
wait
