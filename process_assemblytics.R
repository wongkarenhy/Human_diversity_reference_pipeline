#!/usr/lib64/R/bin/Rscript
## ---------------------------------------------------------------------------------
##
## Script name: process_assemblytics.R
##
## Purpose of script: Cleans up Assemblytics output files
##
## Author: Karen Wong
##
## Date Created: 06/06/2019
##
## Email: karen.wong4@ucsf.edu
##
## ---------------------------------------------------------------------------------
##
## Notes: This script has three parts: filter based on general features, add ngaps and bionano info, filter based on ngaps
##   
##
##    Part 1: Filter based on general features
##      1. segdups.bedpe (defined by 10X Genomics)
##      2. sv_blacklist.bed (defined by 10X Genomics)
##      3. primary chromosome (chr1-22,X,Y)
##      4. ref_gap_size is <-7kb or >1kb
##      5. redundant insertions (INS identified by different scaffolds but same ref position)
##
##    Part 2: Annotate bionano SV calls, adjust query coordinates, and overlap ngaps
##      1. overlap bionano SV calls 
##      2. adjust assm_coords : within script neg strand somehow yielded negative ranges
##      2. adjust query start and end coordinates (only for the INS with negative ref_gap_size)
##        i. remove sequence if start and end coordinates are outside assembly coordinates
##      3. determine whether INS and INS boundaries are overlapping ngaps
##        
##    Part 3: Filter based on ngaps
##      Only INS fulfilling either of the following criteria are kept:
##      1. ngap ==0
##      2. (ngap > 0) & (q_gap_size - ngap >= 50) & (ngap < 10000)
##        
##
## ---------------------------------------------------------------------------------
##

# remove all previous variables
rm(list=ls())

# load up packages  
suppressMessages(library(optparse))
suppressMessages(library(stringr))
suppressMessages(library(GenomicRanges))
suppressMessages(library(data.table))
suppressMessages(library(parallel))

# command line options
## ---------------------------------------------------------------------------------

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

## ---------------------------------------------------------------------------------

# parse user input
threads = opt$threads
dir = opt$dir
BN_path = opt$bn_path
ngap_path = opt$ngap

# load up functions 
source(paste0(dir, "/scripts/insertion_filtering_functions.R"))

# read sample metadata (so we only use samples in the metadata)
metadata = read.table(paste0(dir, "/TMP_sample_metadata.txt"), stringsAsFactors = F)
colnames(metadata) = c("SAMPLE", "SEX", "FASTQ_DIR", "LONGRANGER_DIR", "ASSM_DIR", "BN_DIR", "ENZYME", "SUPERNOVA_VER", "ALT_NAME", "NUCMER_DIR", "POPULATION", "SRC", "PROJECT")

# grab all files in the assemblytics output directory
assemblytics_list = list.files(paste0(dir, "/assemblytics/"), pattern = "filtered_variants.bed")

# only keep samples that are in the metadata sample list
assemblytics_list_sample = str_split_fixed(assemblytics_list, "_|\\.",2)[,1]
assemblytics_list = assemblytics_list[(assemblytics_list_sample %in% metadata$SAMPLE) | (assemblytics_list_sample %in% metadata$ALT_NAME)]

# get a list of samples with BN data
BN_list = list.files(BN_path)
BN_list = str_split_fixed(BN_list, "\\.|_", 2)[,1]

# read the fasta idx file
faidx = fread(paste0(dir,"/discovery/tmp_idx.txt"), stringsAsFactors = F)
colnames(faidx) = c("assm_id", "contig_length", "sample", "haplo")

# read segdups and blacklist files 
segdup = read.table(paste0(dir, "/segdups.bedpe"), stringsAsFactors = F)
sv_bl = read.table(paste0(dir, "/sv_blacklist.bed"), stringsAsFactors = F)

# put them into grange objects
segdup.gr = makeGRangesFromDataFrame(segdup, seqnames.field = "V1", start.field = "V2", end.field = "V3")
sv_bl.gr = makeGRangesFromDataFrame(sv_bl, seqnames.field = "V1", start.field = "V2", end.field = "V3")

# define primary chr_list
chr_list = c(paste0("chr",c(1:22, "X", "Y")))

processAlignment = function(i) {
  
  print(i)
  
  # get the sample name of the file we're handling
  sample = str_split_fixed(i, "_|\\.", 2)[,1]
  
  # get the haplo 
  haplo = substr((str_split_fixed(i, "_", 2)[,2]), 1,3)
  if (haplo!="2.1" & haplo!="2.2"){ #PB samples are unphased
      haplo = "unphased"
  }
  
  # read the assemblytics file
  assemblytics = read.table(paste0(dir, "/assemblytics/",i), stringsAsFactors = F)
  colnames(assemblytics) = c("ref_chr", "ref_start", "ref_end", "assemblytics_id","insert_size", "strand","type", "ref_gap_size" ,"q_gap_size" ,"assm_coords", "method", "adjusted_coords")
  
  
  ## ---------------------------------------------------------------------------------
  ## Part 1: Filter based on general features
  
  
  # remove INS if ref_gap_size is <-7kb or >1kb
  # assemblytics = assemblytics[(assemblytics$ref_gap_size>=(-7000) & assemblytics$ref_gap_size<=1000),]

  # calculate ref_gap_size to query_gap_size ratio
  assemblytics$gap_ratio = assemblytics$ref_gap_size/assemblytics$q_gap_size
  
  # keep only primary chromosomes
  assemblytics = assemblytics[assemblytics$ref_chr %in% chr_list,]
  
  # make a GRange object for assemblytics
  assemblytics.gr = makeGRangesFromDataFrame(assemblytics, seqnames.field = "ref_chr", start.field = "ref_start", end.field = "ref_end")
  
  # filter based on segdup and sv blacklist
  assemblytics = removeSegdupSVblacklist(assemblytics, assemblytics.gr, segdup.gr, sv_bl.gr)
  
  # remove redundant insertions
  assemblytics = removeRedundantInsertion(assemblytics)
  
  
  ## ---------------------------------------------------------------------------------
  ## Part 2: Annotate bionano SV calls, adjust query coordinates, and overlap ngaps
  
  
  # read BN smap if file exists
  if (!(sample %in% BN_list)){ # if no optimal map for this particular sample
    BN_sample = NULL 
  } else {
    BN_sample = readBN(sample, BN_path)
  }
  
  # overlap BN SV size
  assemblytics.gr = makeGRangesFromDataFrame(assemblytics, seqnames.field = "ref_chr", start.field = "ref_start", end.field = "ref_end") # Need a new GRange object 
  assemblytics = overlapBN(BN_sample, assemblytics, assemblytics.gr)
  
  # add BN_validate column based on BN calls
  assemblytics = BN_validate(assemblytics)
  
  # adjust assm_coords for those with negative width
  assemblytics = processAssmCoords(assemblytics)
  
  # adjust query start and end coordinates
  assemblytics = processAdjustedCoords(assemblytics)
  
  # remove sequence if start and end coordinates are outside assembly scaffold coordinates
  assemblytics = checkAssmCoordsWithinLen(assemblytics, faidx, sample, haplo)
  
  # overlap assemblytics INS with assembly ngaps
  if (haplo == "unphased") {
    ngapFileName = paste0(ngap_path, "/",sample, ".ngaps.bed")
    ngapFileSize = file.info(ngapFileName)$size
    
    if (ngapFileSize==0){
      ngapFile = NULL
    } else {
      ngapFile = read.table(ngapFileName, stringsAsFactors = F)
    }
  } else {
    ngapFileName = paste0(ngap_path, "/",sample,"_pseudohap",haplo,".ngaps.bed")
    ngapFile = read.table(ngapFileName, stringsAsFactors = F)
  }
  
  # actual overlap of Ngap
  assemblytics = overlapNgap(assemblytics, ngapFile)

  
  ## ---------------------------------------------------------------------------------
  ## Part 3: Filter based on ngaps
  
  
  # Remove if ngap size is large and has no unique sequence
  assemblytics = assemblytics[(assemblytics$ngap == 0) | (assemblytics$ngap > 0 & (assemblytics$q_gap_size - assemblytics$ngap >= 50) & assemblytics$ngap < 10000) | (assemblytics$ngap_perct<0.5),]

  
  ## ---------------------------------------------------------------------------------
  ## Clean up dataframe for writing output file
  
  
  # Remove unnecessary columns, add sample and haplo info
  assemblytics = assemblytics[order(assemblytics$ref_chr, assemblytics$ref_start),-which(names(assemblytics) %in% c("assemblytics_id", "type", "boundary_left_start", "boundary_left_end", "boundary_right_start", "boundary_right_end"))]
  assemblytics$sample = sample
  assemblytics$haplo = haplo
  
  # Sort and reorder dataframe
  assemblytics = assemblytics[order(assemblytics$ref_chr, assemblytics$ref_start), ] 
  assemblytics = assemblytics[,c("ref_chr","ref_start","ref_end","insert_size","strand","ref_gap_size","q_gap_size","assm_coords","assm_id","assm_start","assm_end","adjusted_coords","adjusted_assm_start","adjusted_assm_end", "adjusted_insert_size", "gap_ratio","BN_size","label_dist","BN_start","BN_end","BN_enzyme","BN_validated","ngap","ngap_boundaries","ngap_boundaries_size_left","ngap_boundaries_size_right", "ngap_perct", "sample","haplo","method")]
  
  # Add INS_id to the last column 
  assemblytics$INS_id = paste0(sample,"_", haplo, "_", 1:nrow(assemblytics))
  
  # Turn off scientific notation
  options(scipen=999)
  write.table(assemblytics, paste0(dir, "/discovery/raw/",sample, "_", haplo,"_assemblytics_raw_results.txt"), col.names = T, row.names = F, quote = F, sep = '\t')
  
  
}

mclapply(assemblytics_list, processAlignment, mc.cores = threads)

#mclapply(assemblytics_list[73], processAlignment, mc.cores = threads)

