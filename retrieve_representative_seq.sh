#!/bin/bash

FILE_TYPE="$1"
WORKDIR="$2"

cd "$WORKDIR"

if [ -f "$WORKDIR"/discovery/final_fasta/"$FILE_TYPE".fa ]; then
    rm "$WORKDIR"/discovery/final_fasta/"$FILE_TYPE".fa
fi

while read -r SAMPLE SEX FASTQ_DIR LONGRANGER_DIR ASSM_DIR BN_DIR ENZYME SUPERNOVA_VER ALT_NAME NUCMER_DIR POPULATION; do
    
    path_sample_name=$(sed 's:.*/::' <<< "$ASSM_DIR" | cut -d_ -f1)

    awk -v SAMPLE="$SAMPLE" '$30==SAMPLE && $5=="-" && $31=="2.1" {print $16":"$23"-"$24}' "$WORKDIR"/discovery/"$FILE_TYPE".txt > "$SAMPLE"_2.1_minus.tmp
    awk -v SAMPLE="$SAMPLE" '$30==SAMPLE && $5=="-" && $31=="2.2" {print $16":"$23"-"$24}' "$WORKDIR"/discovery/"$FILE_TYPE".txt > "$SAMPLE"_2.2_minus.tmp
    awk -v SAMPLE="$SAMPLE" '$30==SAMPLE && $5=="+" && $31=="2.1" {print $16":"$23"-"$24}' "$WORKDIR"/discovery/"$FILE_TYPE".txt > "$SAMPLE"_2.1_plus.tmp
    awk -v SAMPLE="$SAMPLE" '$30==SAMPLE && $5=="+" && $31=="2.2" {print $16":"$23"-"$24}' "$WORKDIR"/discovery/"$FILE_TYPE".txt > "$SAMPLE"_2.2_plus.tmp

    for haplo in 2.1 2.2; do
    
        xargs samtools_1.9 faidx --reverse-complement --mark-strand no "$ASSM_DIR"/"$path_sample_name"_pseudohap"$haplo".fasta < "$SAMPLE"_"$haplo"_minus.tmp >> "$WORKDIR"/discovery/final_fasta/"$FILE_TYPE".fa
        xargs samtools_1.9 faidx "$ASSM_DIR"/"$path_sample_name"_pseudohap"$haplo".fasta < "$SAMPLE"_"$haplo"_plus.tmp >> "$WORKDIR"/discovery/final_fasta/"$FILE_TYPE".fa

    done

    # Remove intermediate files
    rm "$SAMPLE"_2.1_minus.tmp "$SAMPLE"_2.2_minus.tmp "$SAMPLE"_2.1_plus.tmp "$SAMPLE"_2.2_plus.tmp
    
done < "$WORKDIR"/ALL_sample_metadata.txt
