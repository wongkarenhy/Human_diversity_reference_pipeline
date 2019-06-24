#!/bin/bash

CHR="$1"
WORKDIR="$2"
HG38_REF="$3"
CORES="$4"

for i in ${CHR}; do
    
    python2.7 ./scripts/multialign.py \
        -i "${WORKDIR}/discovery/assemblytics_combined_results_with_component_group_big_chr${i}.txt" \
        -r "${HG38_REF}" \
        -t "${WORKDIR}/tmp/" \
        -m "${WORKDIR}/TMP_sample_metadata.txt" \
        -n "${CORES}"

    python2.7 ./scripts/multialign.py \
        -i "${WORKDIR}/discovery/assemblytics_combined_results_with_component_group_small_chr${i}.txt" \
        -r "${HG38_REF}" \
        -t "${WORKDIR}/tmp/" \
        -m "${WORKDIR}/TMP_sample_metadata.txt" \
        -n "${CORES}" \
        -a 0
    
done 
