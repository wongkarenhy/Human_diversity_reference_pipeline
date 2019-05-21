# Need to clean up the data set before we can make comparisons between different technologies
# From the assemblytics output, we need to first remove inserted sequences with ngaps 
#!/usr/lib64/R/bin/Rscript
# Remove all previous variables
rm(list=ls())
######################################################################################################################
#command line options
library(optparse)

option_list = list(
  
  make_option(c("-t", "--threads"), action="store", default=NA, type='numeric',
              help="number of threads"),
  make_option(c("-d", "--dir"), action="store", default=NA, type='character',
              help="work directory"),
  make_option(c("-b", "--bn_path"), action="store", default=NA, type='character',
              help="path to the BN SV7989 folder"),
  make_option(c("-n", "--ngap"), action="store", default=NA, type='character',
              help="path to assembly ngap folder")
  
)
opt = parse_args(OptionParser(option_list=option_list))
################################################################################################################
suppressMessages(library(stringr))
suppressMessages(library(GenomicRanges))

# Parse the user input
threads = opt$threads
dir = opt$dir
BN_path = opt$bn_path
ngap = opt$ngap

source(paste0(dir, "/scripts/insertion_filtering_functions.R"))

# Read sample metadata (so we only use samples in the metadata)
metadata = read.table(paste0(dir, "/TMP_sample_metadata.txt"), stringsAsFactors = F)
colnames(metadata) = c("SAMPLE", "SEX", "FASTQ_DIR", "LONGRANGER_DIR", "ASSM_DIR", "BN_DIR", "ENZYME", "SUPERNOVA_VER", "ALT_NAME", "NUCMER_DIR", "POPULATION", "SRC")

# Grab all files in the assemblytics directory
file_list = list.files(paste0(dir, "/assemblytics/"), pattern = "variants_between_alignments_new.bed")

# Keep samples that are in the metadata sample list
file_list_sample = str_split_fixed(file_list, "_|\\.",2)[,1]
file_list = file_list[(file_list_sample %in% metadata$SAMPLE) | (file_list_sample %in% metadata$ALT_NAME)]

# Get a list of samples with BN DLE data
BN_list = list.files(BN_path)
BN_list = str_split_fixed(BN_list, "\\.", 2)[,1]

# Use the blacklist to remove potentially spurious alignments
segdup = read.table(paste0(dir, "/segdups.bedpe"), stringsAsFactors = F)
sv_bl = read.table(paste0(dir, "/sv_blacklist.bed"), stringsAsFactors = F)

segdup.gr = makeGRangesFromDataFrame(segdup, seqnames.field = "V1", start.field = "V2", end.field = "V3")
sv_bl.gr = makeGRangesFromDataFrame(sv_bl, seqnames.field = "V1", start.field = "V2", end.field = "V3")

chr_list = c(paste0("chr",c(1:22, "X", "Y")))

processAlignment = function(i) {
  
  print(i)
  
  sample = str_split_fixed(i, "_|\\.", 2)[,1]

  haplo = substr((str_split_fixed(i, "_", 2)[,2]), 1,3)
  if (haplo!="2.1" & haplo!="2.2"){ #PB samples are unphased
      haplo = "unphased"
  }

  assemblytics = read.table(paste0(dir, "/assemblytics/",i), stringsAsFactors = F)
  assemblytics = assemblytics[,-11]
  colnames(assemblytics) = c("ref_chr", "ref_start", "ref_end", "assemblytics_id","insert_size", "strand","type", "ref_gap_size" ,"q_gap_size" ,"asm_coords", "adjusted_coords")
  
  # Only keep insertional sequences longer than 50bp
  assemblytics = assemblytics[(assemblytics$type=="Insertion" | assemblytics$type=="Tandem_expansion" | assemblytics$type=="Repeat_expansion") & (assemblytics$insert_size >= 50),]

  # Remove INS if ref_gap_size is <-7kb or >1kb
  assemblytics = assemblytics[(assemblytics$ref_gap_size>=(-7000) & assemblytics$ref_gap_size<=1000),]

  # Calculate ref_gap_size to query_gap_size ratio
  assemblytics$gap_ratio = assemblytics$ref_gap_size/assemblytics$q_gap_size
  
  # Keep only primary chromosomes
  assemblytics = assemblytics[assemblytics$ref_chr %in% chr_list,]
  
  # Make a GRange object for assemblytics
  assemblytics.gr = makeGRangesFromDataFrame(assemblytics, seqnames.field = "ref_chr", start.field = "ref_start", end.field = "ref_end")
  
  # Check segdup 
  assemblytics = removeSegdup(assemblytics, assemblytics.gr, segdup.gr, sv_bl.gr)
  # Read BN smap if it exists
  if (!(sample %in% BN_list)){ # if no optimal map for this particular sample
    BN_sample = NULL 
  } else {
    BN_sample = readBN(sample, BN_path)
  }
  # Check BN SV size
  assemblytics.gr = makeGRangesFromDataFrame(assemblytics, seqnames.field = "ref_chr", start.field = "ref_start", end.field = "ref_end") # Need a new GRange object 
  assemblytics = overlapBN(BN_sample, assemblytics, assemblytics.gr)
  
  # Check ngap 
  # Need to define 2 types of ngap: ngaps_within and ngaps_at_boundaries
  if (haplo == "unphased") {
    ngapFileName = paste0(ngap, "/",sample, ".ngaps.bed")
    ngapFileSize = file.info(ngapFileName)$size
    
    if (ngapFileSize==0){
      ngapFile = NULL
    } else {
      ngapFile = read.table(ngapFileName, stringsAsFactors = F)
    }
  } else {
    ngapFileName = paste0(ngap, "/",sample,"_pseudohap",haplo,".ngaps.bed")
    ngapFile = read.table(ngapFileName, stringsAsFactors = F)
  }
  
  if (!is.null(ngapFile)){
    colnames(ngapFile) = c("ref_assm", "start", "end", "ngap_size")
    ngap.gr = makeGRangesFromDataFrame(ngapFile, seqnames.field = "ref_assm", start.field = "start", end.field = "end")
  }

  # Now we need to overlap assemblytics and ngap
  assemblytics$assm_id = str_split_fixed(assemblytics$asm_coords, ":|-", 4)[,1]
  assemblytics$assm_start = as.numeric(str_split_fixed(assemblytics$asm_coord, ":|-", 4)[,2]) 
  assemblytics$assm_end = as.numeric(str_split_fixed(assemblytics$asm_coord, ":|-", 4)[,3]) 
  assemblytics$boundary_left_start = assemblytics$assm_start - 10
  assemblytics$boundary_left_end = assemblytics$assm_start + 10
  assemblytics$boundary_right_start = assemblytics$assm_end - 10
  assemblytics$boundary_right_end = assemblytics$assm_end + 10
  
  if (is.null(ngapFile)){
    
    assemblytics$ngap = 0
    assemblytics$ngap_boundaries = 0
    assemblytics$ngap_boundaries_size_left = 0
    assemblytics$ngap_boundaries_size_right = 0
    
  } else {
    
    assemblytics_ngap.gr = makeGRangesFromDataFrame(assemblytics, seqnames.field = "assm_id", start.field = "assm_start", end.field = "assm_end")
    assemblytics_boundary_left.gr = makeGRangesFromDataFrame(assemblytics, seqnames.field = "assm_id", start.field = "boundary_left_start", end.field = "boundary_left_end")
    assemblytics_boundary_right.gr = makeGRangesFromDataFrame(assemblytics, seqnames.field = "assm_id", start.field = "boundary_right_start", end.field = "boundary_right_end")
    
    # Actual overlap of Ngap
    assemblytics = overlapNgap(assemblytics, ngapFile, assemblytics_ngap.gr, ngap.gr, assemblytics_boundary_left.gr, assemblytics_boundary_right.gr)
  }
  
  # Remove redundant insertions
  assemblytics = removeRedundantInsertion(assemblytics)
  
  # Remove if ngap size is large and has no unique sequence
  #assemblytics = assemblytics[(assemblytics$ngap == 0) | (assemblytics$ngap > 10 & assemblytics$ngap < 1000 & assemblytics$q_gap_size - assemblytics$ngap >= 150) |(assemblytics$ngap >= 1000 & assemblytics$q_gap_size - assemblytics$ngap >= 1000 & assemblytics$BN_size > 0),]
  assemblytics = assemblytics[(assemblytics$ngap == 0) | (assemblytics$ngap > 0 & (assemblytics$q_gap_size - assemblytics$ngap >= 50) & assemblytics$ngap < 10000),]

  # Remove unnecessary columns, add sample and haplo info
  assemblytics = assemblytics[order(assemblytics$ref_chr, assemblytics$ref_start),-which(names(assemblytics) %in% c("assemblytics_id", "type", "boundary_left_start", "boundary_left_end", "boundary_right_start", "boundary_right_end"))]
  assemblytics$sample = sample
  assemblytics$haplo = haplo
  
  # Parse the adjusted_coordinates column
  # **1. Remove INS if adjusted_coordinates can't be computed given ref_gap_size is negative
  discard = grep("\\[", assemblytics$adjusted_coords)
  if (length(discard)>0){
      assemblytics = assemblytics[-discard,]
  }
  assemblytics$adjusted_assm_start = NA
  assemblytics$adjusted_assm_end = NA

  # **2. Split assemblytics into 2 dataframes (negative ref_gap_size and positive ref_gap_size)
  # we will only use the adjusted_coordinates column for neg_ref_gap
  neg_ref_gap = assemblytics[assemblytics$ref_gap_size<0,]
  pos_ref_gap = assemblytics[assemblytics$ref_gap_size>=0,]

  # **2. Make sure only two sets of coordinates are reported
  # in neg_ref_gap, we need to make sure assm_start and assm_end are found in the adjusted_coordinates column
  neg_ref_gap_list = strsplit(neg_ref_gap$adjusted_coords, ";")
  discard = NULL
  for (s in 1:length(neg_ref_gap_list)){
    
      startPattern = neg_ref_gap$assm_start[s]
      endPattern = neg_ref_gap$assm_end[s]
      
      adjusted_assm_start = str_split_fixed(neg_ref_gap_list[[s]][grep(startPattern,neg_ref_gap_list[[s]])], "-", 2)
      if (nrow(adjusted_assm_start)==1) neg_ref_gap$adjusted_assm_start[s] = adjusted_assm_start[which(adjusted_assm_start != startPattern)]
      else if (nrow(adjusted_assm_start)==2) {
        neg_ref_gap$adjusted_assm_start[s] = min(adjusted_assm_start[which(!adjusted_assm_start %in% startPattern)])
        neg_ref_gap$adjusted_assm_end[s] = max(adjusted_assm_start[which(!adjusted_assm_start %in% startPattern)])
        next
      }
      else discard = c(discard, s)
      
      adjusted_assm_end = str_split_fixed(neg_ref_gap_list[[s]][grep(endPattern,neg_ref_gap_list[[s]])], "-", 2)
      if (nrow(adjusted_assm_end)==1) neg_ref_gap$adjusted_assm_end[s] = adjusted_assm_end[which(adjusted_assm_end != endPattern)]
      else discard = c(discard, s)
  }
  neg_ref_gap = neg_ref_gap[-discard,]

  neg_ref_gap$tmp_s = neg_ref_gap$adjusted_assm_start
  neg_ref_gap$tmp_e = neg_ref_gap$adjusted_assm_end
  
  neg_ref_gap$adjusted_assm_end[(neg_ref_gap$tmp_e < neg_ref_gap$tmp_s)] = neg_ref_gap$tmp_s[(neg_ref_gap$tmp_e < neg_ref_gap$tmp_s)]
  neg_ref_gap$adjusted_assm_start[(neg_ref_gap$tmp_e < neg_ref_gap$tmp_s)] = neg_ref_gap$tmp_e[(neg_ref_gap$tmp_e < neg_ref_gap$tmp_s)]
  neg_ref_gap = subset(neg_ref_gap, select = -c(tmp_s, tmp_e))
  
  # **3. If pos_ref_gap, adjusted_assm_start and adjusted_assm_end are identical to assm_start and assm_end
  pos_ref_gap$adjusted_assm_start = pos_ref_gap$assm_start
  pos_ref_gap$adjusted_assm_end = pos_ref_gap$assm_end

  # **4. Put the two lists together
  assemblytics = rbind.data.frame(pos_ref_gap, neg_ref_gap, stringsAsFactors = F)
  
  # Sort dataframe
  assemblytics = assemblytics[order(assemblytics$ref_chr, assemblytics$ref_start),]
  
  # Add BN_validate column based on BN calls
  assemblytics = BN_validate(assemblytics)
  
  # Add INS_id to the last column 
  assemblytics$INS_id = paste0(sample,"_", haplo, "_", 1:nrow(assemblytics))
  
  # Turn off scientific notation
  options(scipen=999)
  write.table(assemblytics, paste0(dir, "/discovery/raw/",sample, "_", haplo,"_assemblytics_raw_results.txt"), col.names = T, row.names = F, quote = F, sep = '\t')
  
}

mclapply(file_list, processAlignment, mc.cores = threads)
#mclapply(file_list[73], processAlignment, mc.cores = threads)

