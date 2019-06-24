#!/bin/bash
CHR="$1"
WORKDIR="$2"
CORES="$3"

# extract sequences and store pids in array
for i in ${CHR}; do

    # extract sequence
    bash ./scripts/retrieve_representative_seq.sh "singleton_chr${i}" "$WORKDIR" "${WORKDIR}/singleton/TMP_chr${i}" &   
    pids[${i}]=$!
    
done

# wait for all pids to finish
for pid in ${pids[*]}; do
    wait $pid
done

for i in ${CHR}; do   
    
    # repeatmasking
    echo [`date +"%Y-%m-%d %H:%M:%S"`] "   * Repeatmasking singletons from chr"$i" "
    RepeatMasker -qq -species human -xsmall -pa "$CORES" -dir "$WORKDIR"/singleton \
        "$WORKDIR"/discovery/final_fasta/singleton_chr"$i".fa
    
    # calculate non-masked sequence composition
    seqtk comp -u "$WORKDIR"/singleton/singleton_chr"$i".fa.masked > "$WORKDIR"/singleton/seqtk_comp_chr"$i".txt

done
