# Insertion_reference

## Required software:
RepeatMasker version open-4.0.7 (Dfam v2.0) <br>
trf (tandem repeat finder)<br>

## Required R libraries:
stringr <br>
GenomicRanges<br>
optparse <br>
igraph <br>
data.table <br>

## Annotation databases:
REFFLAT file: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz (Version 03312019)<br>
GWAS file: downloaded from the GWAS catalog website but the file is updated weekly. Version used: 02032019

## Other input requirements:
1. Metadata files (Filename has to be **10X_sample_metadata.txt** for 10X data and **PB_sample_metadata.txt** for PB)
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

**Supernova pseudohaplotype unzipped FASTA files** <br>
File naming format: $SAMPLE_pseudohap2.1.fasta AND $SAMPLE_pseudohap2.2.fasta

**Nucmer delta files** <br>
File naming format: OUT_$SAMPLE_pseudohap2.1.delta AND OUT_$SAMPLE_pseudohap2.2.delta

**BioNano SV smap** <br>
Rerun Bionano's run_SV.py for all the samples. <br>
Make sure all smaps are stored in one directory <br>
File naming format: $SAMPLE.smap 

## Usage:<br>
Place reference_config.sh in WORKDIR and change all the variables accordingly <br>
Put segdups.bedpe and sv_blacklist.bed in WORKDIR <br>
Hg38 reference FASTA header has to follow EXACTLY the following format: <br>
**>chrX$** (no space or tab allowed) <br>

**To execute the pipeline, run:** <br>
bash pipeline.sh <br>

**Pipeline input/output directories**<br>
├── 10X_sample_metadata.txt --> **pipeline input**<br>
├── assemblytics --> assemblytics output directory <br>
├── discovery<br>
│   ├── assemblytics_combined_results.txt<br>
│   ├── assemblytics_combined_results_with_component_group.txt<br>
│   ├── assemblytics_component_edge.txt<br>
│   ├── assemblytics_representative_seq_annotated.txt<br>
│   ├── assemblytics_representative_seq_conf_annotated.txt --> input to the make reference pipeline<br>
│   ├── assemblytics_representative_seq_annotated.txt<br>
│   ├── assemblytics_representative_seq_conf_annotated.txt<br>
│   ├── assemblytics_representative_seq.txt<br>
│   ├── comp_repeat.txt<br>
│   ├── final_fasta<br>
│   │   ├── assemblytics_representative_seq.fa<br>
│   │   ├── assemblytics_representative_seq.fa.2.7.7.80.10.50.2000.mask<br>
│   │   ├── repeats<br>
│   │   ├── rep_seq_count.txt<br>
│   │   ├── trf_rep_seq_count.txt<br>
│   │   └── trf_rep_seq.out<br>
│   ├── raw --> output of process_assemblytics.R<br>
│   └── tmp_idx.txt --> output of compile_fasta_idx.sh <br>
├── log<br>
├── multi_results.csv --> output of multialign.py<br>
├── PB_sample_metadata.txt --> **pipeline input**<br>
├── reference_config.sh --> **pipeline input**<br>
├── scripts --> **pipeline input** *all scripts should be stored here* <br>
├── segdups.bedpe --> **pipeline input**<br>
├── sv_blacklist.bed --> **pipeline input**<br>
├── tmp --> output of multialign.py<br>





