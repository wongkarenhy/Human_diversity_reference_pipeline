#!/bin/bash
INPUT="$1"
WORKDIR="$2"
OUTPUT="$3"
CORES="$4"

while read -r LINE; do
    
    SAMPLE=$(awk '{print $28}' <<< "$LINE")
    HAPLO=$(awk '{print $29}' <<< "$LINE")
    ASSM_ID=$(awk '{print $9}' <<< "$LINE")
    ASSM_START=$(awk '{print $13}' <<< "$LINE")
    ASSM_END=$(awk '{print $14}' <<< "$LINE")
    STRAND=$(awk '{print $5}' <<< "$LINE")
    INS_ID=$(awk '{print $32}' <<< "$LINE")
    
    if [ "$HAPLO" = "unphased" ]; then        
        ASSM_DIR=$(grep "$SAMPLE" "$WORKDIR"/TMP_sample_metadata.txt | awk '{print $5}')   
    else
        ASSM_DIR=$(grep "$SAMPLE" "$WORKDIR"/TMP_sample_metadata.txt | awk '{print $5}')
        path_sample_name=$(sed 's:.*/::' <<< "$ASSM_DIR" | cut -d_ -f1)
        ASSM_DIR="$ASSM_DIR"/"$path_sample_name"_pseudohap"$HAPLO".fasta
    fi
    
    if [ "$STRAND" = "-" ]; then
        samtools faidx --reverse-complement --mark-strand no "$ASSM_DIR" "$ASSM_ID":"$ASSM_START"-"$ASSM_END" | tr [a-z] [A-Z] > "$OUTPUT"/"$INS_ID".fa   
    else
        samtools faidx "$ASSM_DIR" "$ASSM_ID":"$ASSM_START"-"$ASSM_END" | tr [a-z] [A-Z] > "$OUTPUT"/"$INS_ID".fa   
    fi
    
    if [ ! -d "$OUTPUT"/"$INS_ID" ]; then
        mkdir "$OUTPUT"/"$INS_ID"
    fi
    
    (RepeatMasker -species human -xsmall -pa 1 -dir "$OUTPUT"/"$INS_ID" "$OUTPUT"/"$INS_ID".fa ;
    trf "$OUTPUT"/"$INS_ID".fa 2 7 7 80 10 50 2000 -m -h ;   
    mv "$INS_ID".fa.2.7.7.80.10.50.2000.mask* "$OUTPUT"/"$INS_ID"/ ;
    seqtk comp "$OUTPUT"/"$INS_ID"/"$INS_ID".fa.2.7.7.80.10.50.2000.mask | awk '{print $1, $2, $9, $9/$2}' > "$OUTPUT"/"$INS_ID"/trf_rep_seq_count.txt ;
    seqtk comp "$OUTPUT"/"$INS_ID".fa | awk '{print $1, ($4+$5)/$2}' > "$OUTPUT"/"$INS_ID"/gc_count.txt ) &
    
    background=( $(jobs -p) )
    if (( ${#background[@]} == CORES )); then
        wait
    fi
  
done < "$INPUT"        
