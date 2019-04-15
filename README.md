# Insertion_reference

# Full pipeline to generate the insertion dataset
**bash /media/KwokRaid05/karen/new_ref/scripts/pipeline.sh CORE WORKDIR** <br>
**For example:** <br>
bash /media/KwokRaid05/karen/new_ref/scripts/pipeline.sh 32 /media/KwokRaid05/karen/new_ref

**Usage:**<br>
Put segdups.bedpe and sv_blacklist.bed in WORKDIR <br>
Open make_annotation.R and change all the database paths under "# Read input sequences and databases"

**Required R libraries:** <br>
stringr <br>
GenomicRanges<br>
optparse <br>
igraph <br>
data.table <br>
