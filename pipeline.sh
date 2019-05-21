#!/bin/bash

# Define variables
CORES="$1"
MODE="$2"
WORKDIR="$3"
ASSEMBLYTICS="$4"
NGAPDIR="$5"
BN_SV7989="$6"
EEE_VCF="$7"

cd "$WORKDIR"

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

if [ ! -d ./discovery/final_fasta ]; then
    mkdir ./discovery/final_fasta
fi
    
if [ ! -d ./discovery/final_fasta/repeats ]; then
    mkdir ./discovery/final_fasta/repeats
fi

# Modify sample metadata depending on the MODE variable
if [ "$MODE" = "ALL" ]; then
    cat "${WORKDIR}/10X_sample_metadata.txt" "${WORKDIR}/PB_sample_metadata.txt" > "${WORKDIR}/TMP_sample_metadata.txt"
elif [ "$MODE" = "10X" ]; then
    cp "${WORKDIR}/10X_sample_metadata.txt" "${WORKDIR}/TMP_sample_metadata.txt"
elif [ "$MODE" = "PB" ]; then
    cp "${WORKDIR}/PB_sample_metadata.txt" "${WORKDIR}/TMP_sample_metadata.txt"
else 
    echo "Error: Incorrect MODE selection"
    exit 1
fi
    
# Calculate ngap size
# This is memory intensive, adjust cores accordingly
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Calculating assembly ngaps"
bash /media/KwokRaid05/karen/new_ref/scripts/calc_ngaps_assembly.sh "$WORKDIR" "$NGAPDIR" 4

# Gather fasta indexes for contig length info
# This outputs "$WORKDIR"/discovery/supernova_idx.txt
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Gathering supernova fasta sequence lengths"
bash ./scripts/compile_fasta_idx.sh "$WORKDIR"
    
# Run assemblytics
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Running assemblytics between alignment"
bash ./scripts/run_assemblytics.sh "$CORES" "$WORKDIR" "$ASSEMBLYTICS" 
    
# Process assemblytics output
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Compiling assemblytics output"
Rscript ./scripts/compile_assemblytics.R -t "$CORES" -d "$WORKDIR" -b "$BN_SV7989" -n "$NGAPDIR"
    
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Combining individual files"
Rscript ./scripts/compile_assemblytics_2.R  -d "$WORKDIR"
    
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Choosing representative insertion sequences"
Rscript ./scripts/find_most_representative_seq_assemblytics.R -d "$WORKDIR" -v "$EEE_VCF" -t "$CORES"

# Extract sequences for TRF
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Extract the insertion sequences"
bash ./scripts/retrieve_representative_seq.sh assemblytics_representative_seq "$WORKDIR"

# Run TRF 
cd ./discovery/final_fasta
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Running tandem repeat finder (TRF)"
trf assemblytics_representative_seq.fa 2 7 7 80 10 50 2000 -f -m -h -d -ngs > trf_rep_seq.out
    
# Ask what percent of the insertions are masked
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Calculating percent masked bases"
seqtk comp assemblytics_representative_seq.fa.2.7.7.80.10.50.2000.mask | \
    awk '{print $1, $2, $9, $9/$2}' > trf_rep_seq_count.txt 
seqtk comp assemblytics_representative_seq.fa | awk '{print $1, $2, $9, $9/$2}' > rep_seq_count.txt

# Do repeatmasking 
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Running repeatMasker"
RepeatMasker --species human -pa 16 -ace -html -dir ./repeats assemblytics_representative_seq.fa

# This script annotates all the raw insertions and filter for a set of high confident SVs
cd "$WORKDIR"
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Annotating insertionns and choosing a high confident set"
Rscript ./scripts/make_annotation.R -d "$WORKDIR"

} # end of pipeline

pipeline 2>&1 | tee $LOGFILE
