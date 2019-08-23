#!/bin/bash                                                                                                                                                                                                 
INPUT="$3"
WORKDIR="$1"
OUTPUT="$2"

while read -r SAMPLE HAPLO ASSM_ID ASSM_START ASSM_END STRAND INS_ID; do

    if [ ! -d "$OUTPUT"/"$INS_ID" ]; then
        mkdir "$OUTPUT"/"$INS_ID"
    elif [ -f "$OUTPUT"/"$INS_ID"/gc_count.txt ] ; then
        continue
    fi

    if [ "$HAPLO" = "unphased" ]; then
        ASSM_DIR=$(grep -w "$SAMPLE" "$WORKDIR"/TMP_sample_metadata.txt | awk '{print $5}')
    else
        ASSM_DIR=$(grep -w "$SAMPLE" "$WORKDIR"/TMP_sample_metadata.txt | awk '{print $5}')
        path_sample_name=$(sed 's:.*/::' <<< "$ASSM_DIR" | cut -d_ -f1)
        ASSM_DIR="$ASSM_DIR"/"$path_sample_name"_pseudohap"$HAPLO".fasta
    fi

    if [ "$STRAND" = "-" ]; then
        samtools faidx --reverse-complement --mark-strand no "$ASSM_DIR" "$ASSM_ID":"$ASSM_START"-"$ASSM_END" | tr [a-z] [A-Z] > "$OUTPUT"/"$INS_ID".fa
    else
        samtools faidx "$ASSM_DIR" "$ASSM_ID":"$ASSM_START"-"$ASSM_END" | tr [a-z] [A-Z] > "$OUTPUT"/"$INS_ID".fa
    fi

    RepeatMasker -species human -xsmall -pa 1 -dir "$OUTPUT"/"$INS_ID" "$OUTPUT"/"$INS_ID".fa ;
    trf "$OUTPUT"/"$INS_ID".fa 2 7 7 80 10 50 2000 -m -h ;
    mv "$INS_ID".fa.2.7.7.80.10.50.2000* "$OUTPUT"/"$INS_ID"/ ;
    seqtk comp "$OUTPUT"/"$INS_ID"/"$INS_ID".fa.2.7.7.80.10.50.2000.mask | awk '{print $1, $2, $9, $9/$2}' > "$OUTPUT"/"$INS_ID"/trf_rep_seq_count.txt ;
    seqtk comp "$OUTPUT"/"$INS_ID".fa | awk '{print $1, ($4+$5)/$2}' > "$OUTPUT"/"$INS_ID"/gc_count.txt 


done < <(tr '!' ' ' <<< $INPUT)

#echo $INPUT
#echo $OUTPUT
