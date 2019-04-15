cores="$1"
workdir="$2"
assemblytics="$3"
metadata="$4"

# Exit immediately upon error
set -e

cd "$workdir"/assemblytics

while read -r SAMPLE SEX FASTQ_DIR LONGRANGER_DIR ASSM_DIR BN_DIR ENZYME SUPERNOVA_VER ALT_NAME NUCMER_DIR POPULATION; do
    
    for haplo in 2.1 2.2; do
    
        # Use Walfred's assemblytics script
        bash "$assemblytics" "$NUCMER_DIR""$haplo".delta \
            "$SAMPLE"_"$haplo" 10000 &
    done
    
    background=( $(jobs -p) )
    if (( ${#background[@]} == cores )); then
        wait
    fi

done < "$metadata"
        
