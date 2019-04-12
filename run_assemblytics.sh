cores="$1"
# Exit immediately upon error
set -e

cd /media/KwokRaid05/karen/new_ref/assemblytics

while read -r SAMPLE SEX FASTQ_DIR LONGRANGER_DIR ASSM_DIR BN_DIR ENZYME SUPERNOVA_VER ALT_NAME NUCMER_DIR POPULATION; do
    
    for haplo in 2.1 2.2; do
    
        # Use Walfred's assemblytics script
        bash /media/KwokRaid02/karen/software/assemblytics_WM/Assemblytics "$NUCMER_DIR""$haplo".delta \
            "$SAMPLE"_"$haplo" 10000 &
    done
    
    background=( $(jobs -p) )
    if (( ${#background[@]} == cores )); then
        wait
    fi

done < /media/KwokRaid05/karen/new_ref/ALL_sample_metadata.txt
        
