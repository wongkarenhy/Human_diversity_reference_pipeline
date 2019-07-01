#!/bin/bash

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

if [ ! -d ./log ]; then
    mkdir ./log
fi

LOGFILE=./log/`date '+%Y-%m-%d'`_pipeline_run.log

# Exit immediately upon error
set -e

## main pipeline
pipeline(){
echo [`date +"%Y-%m-%d %H:%M:%S"`] "#> START: " $0 

if [ ! -d ./assemblytics ]; then
    mkdir ./assemblytics
fi

if [ ! -d ./discovery ]; then
    mkdir ./discovery
fi

if [ ! -d ./discovery/raw ]; then
    mkdir ./discovery/raw
fi

if [ ! -d ./discovery/singleton ]; then
    mkdir ./discovery/singleton
fi

if [ ! -d ./discovery/final_fasta ]; then
    mkdir ./discovery/final_fasta
fi
    
# Modify sample metadata depending on the MODE variable
if [ "$MODE" = "ALL" ]; then
    cp "${WORKDIR}/sample_metadata.txt" "${WORKDIR}/TMP_sample_metadata.txt"
elif [ "$MODE" = "10X" ]; then
    grep -w 10X "${WORKDIR}/sample_metadata.txt" > "${WORKDIR}/TMP_sample_metadata.txt"
elif [ "$MODE" = "PB" ]; then
    grep -w PB "${WORKDIR}/sample_metadata.txt" > "${WORKDIR}/TMP_sample_metadata.txt"
else 
    echo "Error: Incorrect MODE selection"
    exit 1
fi
    
# Calculate ngap size
# This is memory intensive, adjust cores accordingly
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Calculating assembly ngaps"
bash ./scripts/calc_ngaps_assembly.sh "$WORKDIR" "$NGAPDIR" 4

# Gather fasta indexes for contig length info
# This outputs "$WORKDIR"/discovery/supernova_idx.txt
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Gathering supernova fasta sequence lengths"
bash ./scripts/compile_fasta_idx.sh "$WORKDIR"
    
# Run assemblytics
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Running assemblytics between alignment"
bash ./scripts/run_assemblytics.sh "$CORES" "$WORKDIR" "$ASSEMBLYTICS"
    
# Process assemblytics output
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Compiling assemblytics output"
Rscript ./scripts/process_assemblytics.R -t "$CORES" -d "$WORKDIR" -b "$BN_SV7989" -n "$NGAPDIR" -c "$CHR"
    
# Group individual assemblytics files
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Combining individual files"
Rscript ./scripts/combine_individual_assemblytics.R -d "$WORKDIR"
    
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Grouping insertions based on overlapping reference coordinates"
Rscript ./scripts/group_insertions.R -d "$WORKDIR" -t "$CORES"
    
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Repeatmasking singletons"
bash ./scripts/RM_singleton.sh "$CHR" "$WORKDIR" "$CORES"
    
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Running multiple alignment for non-singleton insertions"
bash ./scripts/run_multialign.sh "$CHR" "$WORKDIR" "$HG38_REF" "$CORES"
  
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Processing multiple alignment results for non-singletons"
Rscript ./scripts/process_multi_results.R -d "$WORKDIR"
  
# Extract repeatMasking results 
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Extracting repeatMasker results for non-singletons"
Rscript ./scripts/extract_repeatmasking_for_rep_seq.R -d "$WORKDIR" -t "$CORES"

# Extract sequences for TRF
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Extract the insertion sequences for non-singletons"
bash ./scripts/retrieve_representative_seq.sh assemblytics_representative_seq "$WORKDIR" "$WORKDIR"
    
# Run TRF 
cd ./discovery/final_fasta
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Running tandem repeat finder (TRF)"
trf assemblytics_representative_seq.fa 2 7 7 80 10 50 2000 -f -m -h -d -ngs > trf_rep_seq.out
    
# Ask what percent of the insertions are masked by TRF
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Calculating percent masked bases"
seqtk comp assemblytics_representative_seq.fa.2.7.7.80.10.50.2000.mask | \
    awk '{print $1, $2, $9, $9/$2}' > trf_rep_seq_count.txt 

# Calculate %GC content 
seqtk comp assemblytics_representative_seq.fa | awk '{print $1, ($4+$5)/$2}' > gc_count.txt

# This script annotates all the raw insertions and filter for a set of high confident SVs
cd "$WORKDIR"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Annotating insertionns"
Rscript ./scripts/make_annotation.R -d "$WORKDIR" -t "$CORES" -r "$REFFLAT" -g "$GWAS"
Rscript ./scripts/make_annotation_for_singleton.R -d "$WORKDIR" -t "$CORES" -r "$REFFLAT" -g "$GWAS" -c "$CHR"

# Choose what insertions should be added to the new reference
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Choosing a high confident set of insertions"
Rscript ./scripts/choose_high_conf_ins.R -d "$WORKDIR"
    
} # end of pipeline

pipeline 2>&1 | tee $LOGFILE
