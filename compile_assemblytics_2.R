#!/usr/lib64/R/bin/Rscript
# This script combines the two pseudohaps into a single file per sample
# Remove all previous variables
rm(list=ls())
######################################################################################################################
#command line options
library(optparse)

option_list = list(
  
  make_option(c("-d", "--dir"), action="store", default=NA, type='character',
              help="work directory")

)
opt = parse_args(OptionParser(option_list=option_list))
################################################################################################################
suppressMessages(library(GenomicRanges))

# Prase user input
dir = opt$dir

source(paste0(dir, "/scripts/insertion_filtering_functions.R"))

file_path = paste0(dir, "/discovery/raw/")
file_list = list.files(file_path, pattern = "raw_results.txt")

# Read all the processed assemblytics files
## Right now there are only 254 files (those missing have no BN files)
assemblytics_ALL = NULL
for (i in file_list){
  
  print(i)
  assemblytics = read.table(paste0(file_path,i), stringsAsFactors = F,header = T)
  # For now we're taking out samples if they have no BN data
  # if (any(is.na(assemblytics$BN_size))) next
  
  # make sure there aren't duplicated entries again
  keep = which(!duplicated(assemblytics[,c("ref_chr", "ref_start", "ref_end")]))
  assemblytics = assemblytics[keep,]
  
  assemblytics_ALL = rbind.data.frame(assemblytics_ALL, assemblytics, stringsAsFactors = F)
}

sample_list = unique(assemblytics_ALL$sample)

assemblytics_ALL_pseudohap_combined = combinePseudohaps(sample_list, assemblytics_ALL)

# Turn off scientific notation
options(scipen=999)
write.table(assemblytics_ALL_pseudohap_combined, paste0(dir, "/discovery/assemblytics_pseudohap_combined_results.txt"), col.names = T, row.names = F, quote = F, sep = '\t')






