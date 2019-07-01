#!/bin/bash
CHR="$1"
WORKDIR="$2"
HG38_REF="$3"
CORES="$4"

for i in ${CHR}; do
    
    python2.7 ./scripts/multialign.py \
        -i "${WORKDIR}/discovery/assemblytics_combined_results_with_component_group_big_chr${i}.txt" \
        -r "${HG38_REF}" \
        -t "${WORKDIR}/tmp_chr${i}/" \
        -m "${WORKDIR}/TMP_sample_metadata.txt" \
        -o "${WORKDIR}/discovery/multi_results_chr${i}.txt" \
        -n "${CORES}"

    python2.7 ./scripts/multialign.py \
        -i "${WORKDIR}/discovery/assemblytics_combined_results_with_component_group_small_chr${i}.txt" \
        -r "${HG38_REF}" \
        -t "${WORKDIR}/tmp_chr${i}/" \
        -m "${WORKDIR}/TMP_sample_metadata.txt" \
        -o "${WORKDIR}/discovery/multi_results_chr${i}.txt" \
        -n "${CORES}" \
        -a 0
    
done 
