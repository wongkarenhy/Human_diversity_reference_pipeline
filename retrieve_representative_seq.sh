#!/bin/bash
FILE_TYPE="$1"
WORKDIR="$2"
TMPDIR="$3"

cd "$WORKDIR"

if [ -f "$WORKDIR"/discovery/final_fasta/"$FILE_TYPE".fa ]; then
    rm "$WORKDIR"/discovery/final_fasta/"$FILE_TYPE".fa
fi

if [ ! -d "$TMPDIR" ]; then
    mkdir "$TMPDIR"
fi

while read -r SAMPLE SEX FASTQ_DIR LONGRANGER_DIR ASSM_DIR BN_DIR ENZYME SUPERNOVA_VER ALT_NAME NUCMER_DIR POPULATION SRC PROJECT; do
    

    if [ "$SRC" = "10X" ]; then
        
        awk -v SAMPLE="$SAMPLE" '$28==SAMPLE && $5=="-" && $29=="2.1" {print $9":"$13"-"$14}' "$WORKDIR"/discovery/"$FILE_TYPE".txt > "$TMPDIR"/"$SAMPLE"_2.1_minus.tmp
        awk -v SAMPLE="$SAMPLE" '$28==SAMPLE && $5=="-" && $29=="2.2" {print $9":"$13"-"$14}' "$WORKDIR"/discovery/"$FILE_TYPE".txt > "$TMPDIR"/"$SAMPLE"_2.2_minus.tmp
        awk -v SAMPLE="$SAMPLE" '$28==SAMPLE && $5=="+" && $29=="2.1" {print $9":"$13"-"$14}' "$WORKDIR"/discovery/"$FILE_TYPE".txt > "$TMPDIR"/"$SAMPLE"_2.1_plus.tmp
        awk -v SAMPLE="$SAMPLE" '$28==SAMPLE && $5=="+" && $29=="2.2" {print $9":"$13"-"$14}' "$WORKDIR"/discovery/"$FILE_TYPE".txt > "$TMPDIR"/"$SAMPLE"_2.2_plus.tmp
    
        path_sample_name=$(sed 's:.*/::' <<< "$ASSM_DIR" | cut -d_ -f1)

        for haplo in 2.1 2.2; do

            cat "$TMPDIR"/"$SAMPLE"_"$haplo"_minus.tmp | xargs samtools faidx --reverse-complement --mark-strand no "$ASSM_DIR"/"$path_sample_name"_pseudohap"$haplo".fasta | tr [a-z] [A-Z] >> "$WORKDIR"/discovery/final_fasta/"$FILE_TYPE".fa 
            cat "$TMPDIR"/"$SAMPLE"_"$haplo"_plus.tmp | xargs samtools faidx "$ASSM_DIR"/"$path_sample_name"_pseudohap"$haplo".fasta | tr [a-z] [A-Z] >> "$WORKDIR"/discovery/final_fasta/"$FILE_TYPE".fa

        done

        # Remove intermediate files
        if [ -f "$TMPDIR"/"$SAMPLE"_2.1_minus.tmp ];then
            rm "$TMPDIR"/"$SAMPLE"_2.1_minus.tmp
        fi

        if [ -f "$TMPDIR"/"$SAMPLE"_2.2_minus.tmp ];then
            rm "$TMPDIR"/"$SAMPLE"_2.2_minus.tmp
        fi
        
        if [ -f "$TMPDIR"/"$SAMPLE"_2.1_plus.tmp ];then
            rm "$TMPDIR"/"$SAMPLE"_2.1_plus.tmp
        fi
        
        if [ -f "$TMPDIR"/"$SAMPLE"_2.2_plus.tmp ];then
            rm "$TMPDIR"/"$SAMPLE"_2.2_plus.tmp
        fi

        
    else 
        
        awk -v SAMPLE="$SAMPLE" '$28==SAMPLE && $5=="-" && $29=="unphased" {print $9":"$13"-"$14}' "$WORKDIR"/discovery/"$FILE_TYPE".txt > "$TMPDIR"/"$SAMPLE"_unphased_minus.tmp
        awk -v SAMPLE="$SAMPLE" '$28==SAMPLE && $5=="+" && $29=="unphased" {print $9":"$13"-"$14}' "$WORKDIR"/discovery/"$FILE_TYPE".txt > "$TMPDIR"/"$SAMPLE"_unphased_plus.tmp

        cat "$TMPDIR"/"$SAMPLE"_unphased_minus.tmp | xargs samtools faidx --reverse-complement --mark-strand no "$ASSM_DIR" | tr [a-z] [A-Z] >> "$WORKDIR"/discovery/final_fasta/"$FILE_TYPE".fa
        cat "$TMPDIR"/"$SAMPLE"_unphased_plus.tmp | xargs samtools faidx "$ASSM_DIR" | tr [a-z] [A-Z] >> "$WORKDIR"/discovery/final_fasta/"$FILE_TYPE".fa
        
        # Remove intermediate files
        if [ -f "$TMPDIR"/"$SAMPLE"_unphased_minus.tmp ];then
            rm "$TMPDIR"/"$SAMPLE"_unphased_minus.tmp
        fi
        
        if [ -f "$TMPDIR"/"$SAMPLE"_unphased_plus.tmp ];then
            rm "$TMPDIR"/"$SAMPLE"_unphased_plus.tmp
        fi
        
    fi
    
done < "$WORKDIR"/TMP_sample_metadata.txt
