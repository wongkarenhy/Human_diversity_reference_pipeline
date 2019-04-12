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
  make_option(c("-c", "--CIAPM"), action="store", default=NA, type='character',
              help="CIAPM ngap path"),
  make_option(c("-g", "--GP"), action="store", default=NA, type='character',
              help="1000GP ngap path")
  
)
opt = parse_args(OptionParser(option_list=option_list))
################################################################################################################
suppressMessages(library(stringr))
suppressMessages(library(GenomicRanges))

# Parse the user input
cores = opt$threads
dir = opt$dir
BN_path = opt$bn_path
CIAPM = opt$CIAPM
GP = opt$GP

source(paste0(dir, "/scripts/insertion_filtering_functions.R"))

file_list = list.files(paste0(dir, "/assemblytics/"), pattern = "variants_between_alignments.bed")

# Get a list of samples with BN DLE data
BN_list = list.files(BN_path)
BN_list = substr(BN_list, 1, 7)

# Use the blacklist to remove potentially spurious alignments
segdup = read.table(paste0(dir, "/segdups.bedpe"), stringsAsFactors = F)
sv_bl = read.table(paste0(dir, "/sv_blacklist.bed"), stringsAsFactors = F)

segdup.gr = makeGRangesFromDataFrame(segdup, seqnames.field = "V1", start.field = "V2", end.field = "V3")
sv_bl.gr = makeGRangesFromDataFrame(sv_bl, seqnames.field = "V1", start.field = "V2", end.field = "V3")

chr_list = c(paste0("chr",c(1:22, "X", "Y")))

processAlignment = function(i) {
  
  print(i)
  
  sample = substr(i,1,7)
  haplo = substr(i,9,11)
  
  assemblytics = read.table(paste0(dir, "/assemblytics/",i), stringsAsFactors = F)
  assemblytics = assemblytics[,-11]
  colnames(assemblytics) = c("ref_chr", "ref_start", "ref_end", "assemblytics_id","insert_size", "strand","type", "ref_gap_size" ,"q_gap_size" ,"asm_coords")
  
  # Only keep insertional sequences
  assemblytics = assemblytics[assemblytics$type=="Insertion" | assemblytics$type=="Tandem_expansion" | assemblytics$type=="Repeat_expansion",]
  # Remove based on other alignment features known to be problematic
  assemblytics = assemblytics[assemblytics$q_gap_size >= 10 & assemblytics$insert_size > 50, ]
  
  # Calculate ref_gap_size to query_gap_size ratio
  assemblytics$gap_ratio = assemblytics$ref_gap_size/assemblytics$q_gap_size
  
  # Remove insertions with unreasonable ratio
  # assemblytics = assemblytics[assemblytics$gap_ratio <= 0.7 & assemblytics$gap_ratio >= -3.33, ]
  
  # Keep only primary chromosomes
  assemblytics = assemblytics[assemblytics$ref_chr %in% chr_list,]
  
  # Make a GRange object for assemblytics
  assemblytics.gr = makeGRangesFromDataFrame(assemblytics, seqnames.field = "ref_chr", start.field = "ref_start", end.field = "ref_end")
  
  # Check segdup 
  assemblytics = removeSegdup(assemblytics, assemblytics.gr, segdup.gr, sv_bl.gr)
  # Read BN smap if it exists
  if (!(sample %in% BN_list)){
    BN_sample = NULL 
  } else {
    BN_sample = readBN(sample, BN_path)
  }
  # Check BN SV size
  assemblytics.gr = makeGRangesFromDataFrame(assemblytics, seqnames.field = "ref_chr", start.field = "ref_start", end.field = "ref_end") # Need a new GRange object 
  assemblytics = overlapBN(BN_sample, assemblytics, assemblytics.gr)
  
  # Check ngap 
  # Need to define 2 types of ngaps_within and ngaps_at_boundaries
  # Read the ngap files (be careful since the ngap files are located at two different paths)
  if (grepl("BC", sample)){ # if CIAPM samples
    ngap = read.table(paste0(CIAPM, "/",sample,"_pseudohap",haplo,".ngaps.bed"), stringsAsFactors = F)
  } else {
    ngap = read.table(paste0(GP ,"/",sample,"_pseudohap",haplo,".ngaps.bed"), stringsAsFactors = F)
  }
  colnames(ngap) = c("ref_assm", "start", "end", "ngap_size")
  ngap.gr = makeGRangesFromDataFrame(ngap, seqnames.field = "ref_assm", start.field = "start", end.field = "end")
  
  # Now we need to overlap assemblytics and ngap
  assemblytics$assm_id = str_split_fixed(assemblytics$asm_coords, ":|-", 4)[,1]
  assemblytics$assm_start = as.numeric(str_split_fixed(assemblytics$asm_coord, ":|-", 4)[,2]) 
  assemblytics$assm_end = as.numeric(str_split_fixed(assemblytics$asm_coord, ":|-", 4)[,3]) 
  assemblytics$boundary_left_start = assemblytics$assm_start - 10
  assemblytics$boundary_left_end = assemblytics$assm_start + 10
  assemblytics$boundary_right_start = assemblytics$assm_end - 10
  assemblytics$boundary_right_end = assemblytics$assm_end + 10
  
  #print("Making ngap GRange objects")
  
  assemblytics_ngap.gr = makeGRangesFromDataFrame(assemblytics, seqnames.field = "assm_id", start.field = "assm_start", end.field = "assm_end")
  assemblytics_boundary_left.gr = makeGRangesFromDataFrame(assemblytics, seqnames.field = "assm_id", start.field = "boundary_left_start", end.field = "boundary_left_end")
  assemblytics_boundary_right.gr = makeGRangesFromDataFrame(assemblytics, seqnames.field = "assm_id", start.field = "boundary_right_start", end.field = "boundary_right_end")
  
  #print("Overlapping Ngap file")
  
  assemblytics = overlapNgap(assemblytics, ngap, assemblytics_ngap.gr, ngap.gr, assemblytics_boundary_left.gr, assemblytics_boundary_right.gr)
  
  #print("Finished overlapping Ngap file")
  
  # Remove redundant insertions
  assemblytics = removeRedundantInsertion(assemblytics)
  
  # Remove if ngap size is large and has no unique sequence
  #assemblytics = assemblytics[(assemblytics$ngap == 0) | (assemblytics$ngap > 10 & assemblytics$ngap < 1000 & assemblytics$q_gap_size - assemblytics$ngap >= 150) |(assemblytics$ngap >= 1000 & assemblytics$q_gap_size - assemblytics$ngap >= 1000 & assemblytics$BN_size > 0),]
  assemblytics = assemblytics[(assemblytics$ngap == 0) | (assemblytics$ngap > 10 & assemblytics$q_gap_size - assemblytics$ngap >= 150),]
  
  # Remove unnecessary columns, add sample and haplo info
  assemblytics = assemblytics[order(assemblytics$ref_chr, assemblytics$ref_start),-which(names(assemblytics) %in% c("assemblytics_id", "type", "boundary_left_start", "boundary_left_end", "boundary_right_start", "boundary_right_end"))]
  assemblytics$sample = sample
  assemblytics$haplo = haplo
  
  write.table(assemblytics, paste0(dir, "/discovery/raw/",sample, "_", haplo,"_assemblytics_raw_results.txt"), col.names = T, row.names = F, quote = F, sep = '\t')
  
}

mclapply(file_list, processAlignment, mc.cores = cores)
#mclapply(file_list[55], processAlignment, mc.cores = cores)
