#!/bin/bash
FILE_TYPE="$1"
WORKDIR="$2"

cd "$WORKDIR"

if [ -f "$WORKDIR"/discovery/final_fasta/"$FILE_TYPE".fa ]; then
    rm "$WORKDIR"/discovery/final_fasta/"$FILE_TYPE".fa
fi

while read -r SAMPLE SEX FASTQ_DIR LONGRANGER_DIR ASSM_DIR BN_DIR ENZYME SUPERNOVA_VER ALT_NAME NUCMER_DIR POPULATION SRC; do
    

    if [ "$SRC" = "10X" ]; then
        
        awk -v SAMPLE="$SAMPLE" '$26==SAMPLE && $5=="-" && $27=="2.1" {print $9":"$13"-"$14}' "$WORKDIR"/discovery/"$FILE_TYPE".txt > "$SAMPLE"_2.1_minus.tmp
        awk -v SAMPLE="$SAMPLE" '$26==SAMPLE && $5=="-" && $27=="2.2" {print $9":"$13"-"$14}' "$WORKDIR"/discovery/"$FILE_TYPE".txt > "$SAMPLE"_2.2_minus.tmp
        awk -v SAMPLE="$SAMPLE" '$26==SAMPLE && $5=="+" && $27=="2.1" {print $9":"$13"-"$14}' "$WORKDIR"/discovery/"$FILE_TYPE".txt > "$SAMPLE"_2.1_plus.tmp
        awk -v SAMPLE="$SAMPLE" '$26==SAMPLE && $5=="+" && $27=="2.2" {print $9":"$13"-"$14}' "$WORKDIR"/discovery/"$FILE_TYPE".txt > "$SAMPLE"_2.2_plus.tmp

    
        path_sample_name=$(sed 's:.*/::' <<< "$ASSM_DIR" | cut -d_ -f1)

        for haplo in 2.1 2.2; do

            xargs samtools faidx --reverse-complement --mark-strand no "$ASSM_DIR"/"$path_sample_name"_pseudohap"$haplo".fasta < "$SAMPLE"_"$haplo"_minus.tmp >> "$WORKDIR"/discovery/final_fasta/"$FILE_TYPE".fa
            xargs samtools faidx "$ASSM_DIR"/"$path_sample_name"_pseudohap"$haplo".fasta < "$SAMPLE"_"$haplo"_plus.tmp >> "$WORKDIR"/discovery/final_fasta/"$FILE_TYPE".fa

        done

        # Remove intermediate files
        rm "$SAMPLE"_2.1_minus.tmp "$SAMPLE"_2.2_minus.tmp "$SAMPLE"_2.1_plus.tmp "$SAMPLE"_2.2_plus.tmp

    else 
        
        awk -v SAMPLE="$SAMPLE" '$26==SAMPLE && $5=="-" && $27=="unphased" {print $9":"$13"-"$14}' "$WORKDIR"/discovery/"$FILE_TYPE".txt > "$SAMPLE"_unphased_minus.tmp
        awk -v SAMPLE="$SAMPLE" '$26==SAMPLE && $5=="+" && $27=="unphased" {print $9":"$13"-"$14}' "$WORKDIR"/discovery/"$FILE_TYPE".txt > "$SAMPLE"_unphased_plus.tmp

        xargs samtools faidx --reverse-complement --mark-strand no "$ASSM_DIR" < "$SAMPLE"_unphased_minus.tmp >> "$WORKDIR"/discovery/final_fasta/"$FILE_TYPE".fa
        xargs samtools faidx "$ASSM_DIR" < "$SAMPLE"_unphased_plus.tmp >> "$WORKDIR"/discovery/final_fasta/"$FILE_TYPE".fa
        
        # Remove intermediate files
        rm "$SAMPLE"_unphased_minus.tmp "$SAMPLE"_unphased_plus.tmp

    fi
    
done < "$WORKDIR"/TMP_sample_metadata.txt
