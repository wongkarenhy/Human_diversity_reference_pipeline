# Insertion_reference

**Required software:**<br>
RepeatMasker version open-4.0.7 (Dfam v2.0) <br>
trf (tandem repeat finder)<br>

**Required R libraries:** <br>
stringr <br>
GenomicRanges<br>
optparse <br>
igraph <br>
data.table <br>

**Annotation databases:**<br>
Gene file: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz (Version 03312019)<br>

Exon file: obtained from UCSC Table Browser (ncbiRefSeq, limited to Exons plus0) 

Non_coding: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/lincRNAsTranscripts.txt.gz (Version 12212015)<br>

gencodeV29_nrRNA file: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.long_noncoding_RNAs.gtf.gz (Version 09292018)<br>
> gunzip gencode.v29.long_noncoding_RNAs.gtf.gz<br>
> cut -f1,4,5 gencode.v29.long_noncoding_RNAs.gtf | grep -v ^# > gencode.v29.long_noncoding_RNAs_simplified.gtf<br>

gencodeV29_pseudo file: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.2wayconspseudos.gtf.gz (Version 09282018)
> gunzip gencode.v29.2wayconspseudos.gtf.gz <br>
> cut -f1,4,5 gencode.v29.2wayconspseudos.gtf | grep -v ^# > gencode.v29.2wayconspseudos_simplified.gtf 

GWAS file: downloaded from the GWAS catalog website but the file is updated weekly. Version used: 02032019

**Major output files**<br>
1. assemblytics_representative_seq_annotated.txt<br>
2. assemblytics_representative_seq_conf_annotated.txt <br>


# Usage:<br>
bash pipeline.sh CORES MODE WORKDIR ASSEMBLYTICS_DIR NGAPDIR BN_SV7989_DIR EEE_VCF_DIR <br>
**Example command:** <br>
bash /media/KwokRaid05/karen/new_ref4/scripts/pipeline.sh 32 10X \ <br>
    /media/KwokRaid05/karen/new_ref4 \ <br>
    /media/KwokRaid02/karen/software/assemblytics_WM/ \ <br>
    /media/KwokRaid04/assembly_ngap/ \ <br>
    /media/KwokRaid02/karen/database/BN_SV7989/ \ <br>
    /media/KwokRaid05/karen/new_ref/published_genomes/EEE/EEE_SV-Pop_1.ALL.sites.20181204.vcf <br>


MODE can be 10X, PB, or ALL <br> 
Put segdups.bedpe and sv_blacklist.bed in WORKDIR <br>
Open make_annotation.R and change all the database paths under "# Read input sequences and databases"
Old reference FASTA header has to follow EXACTLY the following format: <br>
**>chrX$** (no space or tab allowed) <br>

**Input to part 1:**<br>
1. Metadata file (Filename has to be **10X_sample_metadata.txt** for 10X data and **PB_sample_metadata.txt** for PB)
2. Supernova pseudohaplotype unzipped FASTA files (2 files per sample if 10X data)
3. Nucmer delta files (2 files per sample)
4. Ngap bed file (2 files per sample)
5. BioNano SV smap (1 file per sample)

# Part 0.1
## Make a sample metadata file
* The metadata file should contain the following columns (* not actually used by pipeline):<br>
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

# Part 0.2
## Run supernova to generate de novo assemblies (FASTA format) <br>
File naming format: $SAMPLE_pseudohap2.1.fasta AND $SAMPLE_pseudohap2.2.fasta

# Part 0.3
## Generate nucmer delta files for all pseudohaplotype assemblies<br>
File naming format: OUT_$SAMPLE_pseudohap2.1.delta AND OUT_$SAMPLE_pseudohap2.2.delta

# Part 0.4
## Generate ngap bed files for all pseudohaplotype assemblies <br>
File naming format: $SAMPLE_pseudohap2.1.ngaps.bed AND $SAMPLE_pseudohap2.2.ngaps.bed

# Part 0.5
## Generate BioNano SV smaps
Rerun Bionano's run_SV.py for all the samples. <br>
All smaps are stored in one directory <br>
File naming format: $SAMPLE.smap 

# Part 1
## Full pipeline to generate the insertion dataset
**bash /media/KwokRaid05/karen/new_ref/scripts/pipeline.sh CORE WORKDIR** <br>
**For example:** <br>
bash /media/KwokRaid05/karen/new_ref/scripts/pipeline.sh 32 /media/KwokRaid05/karen/new_ref

# Part 2
## Construct new reference
> bash /path/to/make_reference.sh WORKDIR NEW_REF_NAME OLD_REF <br>

**For example:** <br>
> bash /media/KwokRaid05/karen/new_ref/scripts/make_reference.sh /media/KwokRaid05/karen/new_ref final_10x \ <br>
    /media/KwokRaid02/karen/reference/hg38/hg38_primary.fa<br>




