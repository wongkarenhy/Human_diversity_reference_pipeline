# Define variables
CORES="$1"
WORKDIR="$2"
ASSEMBLYTICS="$3"
METADATA="$4"

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

# Gather fasta indexes for contig length info
# This outputs "$WORKDIR"/discovery/supernova_idx.txt
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Gathering supernova fasta indexes"
bash /media/KwokRaid05/karen/new_ref/scripts/compile_fasta_idx.sh "$WORKDIR" "$METADATA"    

# Run assemblytics
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Running assemblytics between alignment"
bash ./scripts/run_assemblytics.sh "$CORES" "$WORKDIR" "$ASSEMBLYTICS" "$METADATA"

# Process assemblytics output
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Compiling assemblytics output"
Rscript ./scripts/compile_assemblytics.R -t "$CORES" -d "$WORKDIR" \
    -b "/media/KwokRaid02/karen/database/BN_SV7989/" \
    -c /media/KwokRaid04/CIAPM/CIAPM_supernova/ngap/ \
    -g /media/KwokRaid04/1000GP/supernova2/1000GP_sn2/ngap/
    
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Combining the two pseudohaplotypes"
Rscript ./scripts/compile_assemblytics_2.R  -d "$WORKDIR" # combine pseudohaplotypes
    
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Choosing representative insertion sequences"
Rscript ./scripts/find_most_representative_seq_assemblytics.R -d "$WORKDIR" \
    -v /media/KwokRaid05/karen/new_ref/published_genomes/EEE/EEE_SV-Pop_1.ALL.sites.20181204.vcf

# Extract sequences for TRF
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Extract the insertion sequences"
bash ./scripts/retrieve_representative_seq.sh \
    assemblytics_representative_seq "$WORKDIR"/discovery

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
echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Annotating insertionns and choosing a high confident set"
Rscript ./scripts/make_annotation.R -d "$WORKDIR"

} # end of pipeline

pipeline 2>&1 | tee $LOGFILE
