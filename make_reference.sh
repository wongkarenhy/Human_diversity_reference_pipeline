#!/bin/bash

# Exit immediately upon error
set -e

usage() {
    NAME=$(basename $0)
    cat <<EOF
EOF
Usage:
  ${NAME} 
Please edit reference_config.sh accordingly before each run. The config file, 
along with 10X_sample_metadata.txt, PB_sample_metadata.txt, segdups.bedpe, 
and sv_blacklist.bed must be placed together with this script in WORKDIR.

EOF
}

## load variables
if [[ ! -f ./reference_config.sh ]]; then
 usage
 echo "Error: Missing configuration file (reference_config.sh) !"
 exit 1
fi

cd "$WORKDIR"

source ./reference_config.sh

# Create reference path
if [ ! -d ./discovery/final_fasta/ref/ ]; then
    mkdir ./discovery/final_fasta/ref/
fi

LOGFILE=log/`date '+%Y-%m-%d'`_make_reference.log

## main pipeline
pipeline(){
echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> START: " $0 

# Retrieve the confident sequence
bash /media/KwokRaid05/karen/new_ref/scripts/retrieve_representative_seq.sh \
    assemblytics_representative_seq_conf_annotated "$WORKDIR"
    
mv ./discovery/final_fasta/assemblytics_representative_seq_conf_annotated.fa ./discovery/final_fasta/ref/"$NEW_REF_NAME".fa
    
# Get the key
awk '{print $1":"$2"-"$3, $9":"$13"-"$14}' ./discovery/assemblytics_representative_seq_conf_annotated.txt > ./discovery/final_fasta/ref/ID_key.txt 

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
    -r "$HG38_REF" \
    -f "$WORKDIR"/discovery/final_fasta/ref/"$NEW_REF_NAME"_name_changed.fa \
    -s hg38_"$NEW_REF_NAME".fa -n record_"$NEW_REF_NAME".tsv > log.txt 

# Index the new reference
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Indexing the new reference file"    
bwa index hg38_"$NEW_REF_NAME".fa

} # end of pipeline

pipeline 2>&1 | tee $LOGFILE

