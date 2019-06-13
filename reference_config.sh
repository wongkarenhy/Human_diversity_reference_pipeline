#---------------------------------------- reference_config.sh -------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------
# Configuration file for the new reference pipeline
# Place this script together with pipeline.sh 

# number of CPUs to be used (not every step will take advantage of multi-threading)
CORES=32

# mode can be ALL, PB, or 10X
MODE='10X'

# all output files will be writting to the workdir
WORKDIR="/media/KwokRaid05/karen/new_ref2/"

# path to the modified assemblytics tool
ASSEMBLYTICS="/media/KwokRaid02/karen/software/assemblytics_WM/"

# path to the ngap files
NGAPDIR="/media/KwokRaid04/assembly_ngap/"

# path to the BN SV files
BN_SV7989="/media/KwokRaid02/karen/database/BN_SV7989/"

# path to the hg38 primary reference
HG38_REF="/media/KwokRaid02/karen/reference/hg38/hg38_primary.fa"

# gene refFlat file for insertion annotations
REFFLAT="/media/KwokRaid02/karen/database/genome_annotation/03312019/refFlat.txt"

# gwas variant list for insertion annotation
GWAS="/media/KwokRaid02/karen/database/gwasCatalog.txt"

# give the reference a name
NEW_REF_NAME="NEW_06112019"
