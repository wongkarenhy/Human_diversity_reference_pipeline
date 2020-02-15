# Human Diversity Reference Pipeline
## By Karen Wong and Walfred Ma

## Required software:
RepeatMasker version open-4.0.7 (Dfam v2.0) <br>
trf (tandem repeat finder) <br>
seqtk v1.0-r82-dirty <br>
kalign 2.0.4 <br>
samtools v1.9 <br>

## Required R libraries:
stringr <br>
GenomicRanges<br>
optparse <br>
igraph <br>
data.table <br>
parallel <br>
plyr <br>

## Other input requirements:
1. Metadata file (Filename has to be **sample_metadata.txt**)
2. Supernova pseudohaplotype unzipped FASTA files (2 files per sample if 10X data)
3. Nucmer delta files (2 files per sample)
4. BioNano SV smap (1 file per sample)

**Metadata file** <br>
The metadata file should contain the following columns (* *not actually used by pipeline*):<br>
1. SAMPLE 
2. SEX *
3. FASTQ_DIR *
4. LONGRANGER_DIR 
5. SUPERNOVA_DIR 
6. BIONANO_ASSEMBLY_DIR *
7. BN_ENZYME *
8. SUPERNOVA_VERSION *
9. ALT_SAMPLE_NAME *
10. NUCMER_DIR 
11. POPULATION 
12. SRC (This can either be 10X or PB) 
13. PROJECT * (samples might have come from different projects)

**Supernova pseudohaplotype unzipped FASTA files** <br>
File naming format: $SAMPLE_pseudohap2.1.fasta AND $SAMPLE_pseudohap2.2.fasta

**Nucmer delta files** <br>
File naming format: OUT_$SAMPLE_pseudohap2.1.delta AND OUT_$SAMPLE_pseudohap2.2.delta

**BioNano SV smap** <br>
Rerun Bionano's run_SV.py for all the samples. <br>
Make sure all smaps are stored in one directory <br>
File naming format: $SAMPLE.smap <br>
Append BN enzyme name to the end of every line in smap <br>

## Usage:<br>
Place reference_config.sh in WORKDIR and change all the variables accordingly <br>
Put segdups.bedpe and sv_blacklist.bed in WORKDIR <br>
Hg38 reference FASTA header has to follow EXACTLY the following format: <br>
**>chrX$** (no space or tab allowed) <br>

**To execute the pipeline, run:** <br>
bash pipeline.sh <br>

**To make the reference, run:** <br>
bash make_reference.sh <br>

