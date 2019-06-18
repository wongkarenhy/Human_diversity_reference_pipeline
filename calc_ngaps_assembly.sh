#!/bin/bash

# Exit immediately upon error
set -e

workdir="$1"
ngapdir="$2"
cores="$3"

while read -r SAMPLE SEX FASTQ_DIR LONGRANGER_DIR ASSM_DIR BN_DIR ENZYME SUPERNOVA_VER ALT_NAME NUCMER_DIR POPULATION SRC PROJECT; do

    if [ "$SRC" = "10X" ]; then 

        # Need to take care of the GM/NA difference
        path_sample_name=$(sed 's:.*/::' <<< "$ASSM_DIR" | cut -d_ -f1)

        if [ -z "$path_sample_name" ]; then
            path_sample_name="$SAMPLE"
        fi

        for i in 2.1 2.2; do
             
            if  [ ! -s "${ngapdir}/${SAMPLE}_pseudohap${i}.ngaps.bed" ]; then

                perl -pe '/^>/ ? print "\n" : chomp' "$ASSM_DIR"/"$path_sample_name"_pseudohap"$i".fasta | \
                    perl -ne 'choqmp;if( />(.*)/){$head = $1; $i=0; next};@a=split("",$_); foreach(@a){$i++;if($_ eq "N" && $s ==0 ){print "$head\t$i"; $s =1}elsif($s==1 && $_ ne "N"){print "\t$i\n";$s=0}}' - | \
                    awk '{if ($8-$7>=1) print $1"\t"$7"\t"$8"\t"$8-$7}' > "${ngapdir}/${SAMPLE}_pseudohap${i}.ngaps.bed" & 
            fi 
            
        done
        
    elif  [ "$SRC" = "PB" ]; then 
        
        if [ ! -f "${ngapdir}/${SAMPLE}.ngaps.bed" ]; then
        
            perl -pe '/^>/ ? print "\n" : chomp' "$ASSM_DIR" | \
                perl -ne 'choqmp;if( />(.*)/){$head = $1; $i=0; next};@a=split("",$_); foreach(@a){$i++;if($_ eq "N" && $s ==0 ){print "$head\t$i"; $s =1}elsif($s==1 && $_ ne "N"){print "\t$i\n";$s=0}}' - | \
                awk '{print $1"\t"$(NF-1)"\t"$NF"\t"$NF-$(NF-1)}' > "${ngapdir}/${SAMPLE}.ngaps.bed" & 
        
        fi
        
    else
        echo "Error: unknwon source!"
        exit 1
    fi
    
    background=( $(jobs -p) )
    if (( ${#background[@]} == cores )); then
        wait
    fi

done < "${workdir}/TMP_sample_metadata.txt"
    
exit 0    
