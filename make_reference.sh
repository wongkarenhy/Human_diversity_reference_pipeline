#!/bin/bash

# Exit immediately upon error
set -e

# Parse user input
WORKDIR="$1"
NEW_REF_NAME="$2"
OLD_REF="$3"

cd "$WORKDIR"

# Create reference path
if [ ! -d ./discovery/final_fasta/ref/ ];then
    mkdir ./discovery/final_fasta/ref/
fi

LOGFILE=log/`date '+%Y-%m-%d'`_make_reference.log

## main pipeline
pipeline(){
echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> START: " $0 

# Retrieve the confident sequence
bash /media/KwokRaid05/karen/new_ref/scripts/retrieve_representative_seq.sh \
    assemblytics_representative_seq_conf2_annotated "$WORKDIR"
    
# if [ -f "$WORKDIR"/discovery/final_fasta/ref/"$NEW_REF_NAME".fa ]; then
#     rm "$WORKDIR"/discovery/final_fasta/ref/"$NEW_REF_NAME".fa
# fi

# # Extract the relevant FASTA sequence
# echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Extracting insertion sequences"    
# while read -r SAMPLE SEX FASTQ_DIR LONGRANGER_DIR ASSM_DIR BN_DIR ENZYME SUPERNOVA_VER ALT_NAME NUCMER_DIR POPULATION; do
    
#     path_sample_name=$(sed 's:.*/::' <<< "$ASSM_DIR" | cut -d_ -f1)

#     awk -v SAMPLE="$SAMPLE" '$22==SAMPLE && $5=="-" && $23!="2.2" && ($58 != "(G)n" || $18 == 0) && ($58 != "(C)n" || $18 == 0) {print $15":"$43"-"$44}' ./assemblytics_representative_seq_conf_annotated.txt > "$SAMPLE"_2.1_minus.tmp
#     awk -v SAMPLE="$SAMPLE" '$22==SAMPLE && $5=="-" && $23=="2.2" && ($58 != "(G)n" || $18 == 0) && ($58 != "(C)n" || $18 == 0) {print $15":"$43"-"$44}' ./assemblytics_representative_seq_conf_annotated.txt > "$SAMPLE"_2.2_minus.tmp
#     awk -v SAMPLE="$SAMPLE" '$22==SAMPLE && $5=="+" && $23!="2.2" && ($58 != "(G)n" || $18 == 0) && ($58 != "(C)n" || $18 == 0) {print $15":"$43"-"$44}' ./assemblytics_representative_seq_conf_annotated.txt > "$SAMPLE"_2.1_plus.tmp
#     awk -v SAMPLE="$SAMPLE" '$22==SAMPLE && $5=="+" && $23=="2.2" && ($58 != "(G)n" || $18 == 0) && ($58 != "(C)n" || $18 == 0) {print $15":"$43"-"$44}' ./assemblytics_representative_seq_conf_annotated.txt > "$SAMPLE"_2.2_plus.tmp

#     for haplo in 2.1 2.2; do
    
#         xargs samtools_1.9 faidx --reverse-complement --mark-strand no "$ASSM_DIR"/"$path_sample_name"_pseudohap"$haplo".fasta < "$SAMPLE"_"$haplo"_minus.tmp >> ./final_fasta/ref/"$NEW_REF_NAME".fa
#         xargs samtools_1.9 faidx "$ASSM_DIR"/"$path_sample_name"_pseudohap"$haplo".fasta < "$SAMPLE"_"$haplo"_plus.tmp >> ./final_fasta/ref/"$NEW_REF_NAME".fa

#     done

#     # Remove intermediate files
#     rm "$SAMPLE"_2.1_minus.tmp "$SAMPLE"_2.2_minus.tmp "$SAMPLE"_2.1_plus.tmp "$SAMPLE"_2.2_plus.tmp
    
# done < /media/KwokRaid05/karen/new_ref/ALL_sample_metadata.txt

mv ./discovery/final_fasta/assemblytics_representative_seq_conf2_annotated.fa ./discovery/final_fasta/ref/"$NEW_REF_NAME".fa
    
# Get the key
awk '{print $1":"$2"-"$3, $16":"$23"-"$24}' ./discovery/assemblytics_representative_seq_conf2_annotated.txt > ./discovery/final_fasta/ref/ID_key.txt 

cd ./discovery/final_fasta/ref/

cp "$NEW_REF_NAME".fa "$NEW_REF_NAME"_name_changed.fa
grep ">" "$NEW_REF_NAME"_name_changed.fa | cut -c2- | grep -Fwf - ID_key.txt > ID_key_subsetted.txt 

# Change the FASTA header
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Changing the FASTA file headers"    
while read -r line; do

    replace=$(awk '{print $2}' <<< "$line")
    substitute=$(awk '{print $1}' <<< "$line")

    sed -i "s/"$replace"/"$substitute"/g" "$NEW_REF_NAME"_name_changed.fa 
    
done < ID_key_subsetted.txt 

echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Making the new reference file"    
python2.7 "$WORKDIR"/scripts/reference_liftover.py \
    -r "$OLD_REF" \
    -f "$WORKDIR"/discovery/final_fasta/ref/"$NEW_REF_NAME"_name_changed.fa \
    -s hg38_"$NEW_REF_NAME".fa -n record_"$NEW_REF_NAME".tsv > log.txt 

# Index the new reference
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Indexing the new reference file"    
bwa index hg38_"$NEW_REF_NAME".fa

} # end of pipeline

pipeline 2>&1 | tee $LOGFILE

