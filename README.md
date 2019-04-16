# Insertion_reference

# Full pipeline to generate the insertion dataset
**bash /media/KwokRaid05/karen/new_ref/scripts/pipeline.sh CORE WORKDIR** <br>
**For example:** <br>
bash /media/KwokRaid05/karen/new_ref/scripts/pipeline.sh 32 /media/KwokRaid05/karen/new_ref

**Usage:**<br>
Put segdups.bedpe and sv_blacklist.bed in WORKDIR <br>
Open make_annotation.R and change all the database paths under "# Read input sequences and databases"

**Software used:**<br>
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
