#!/usr/lib64/R/bin/Rscript
# This script combines all individual processed assemblytics files
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
suppressMessages(library(stringr))

# Prase user input
dir = opt$dir

source(paste0(dir, "/scripts/insertion_filtering_functions.R"))

metadata = read.table(paste0(dir, "/TMP_sample_metadata.txt"), stringsAsFactors = F)
colnames(metadata) = c("SAMPLE", "SEX", "FASTQ_DIR", "LONGRANGER_DIR", "ASSM_DIR", "BN_DIR", "ENZYME", "SUPERNOVA_VER", "ALT_NAME", "NUCMER_DIR", "POPULATION", "SRC")

file_path = paste0(dir, "/discovery/raw/")
file_list = list.files(file_path, pattern = "raw_results.txt")

# Keep samples that are in the metadata sample list
file_list_sample = str_split_fixed(file_list, "_",2)[,1]
file_list = file_list[(file_list_sample %in% metadata$SAMPLE) | (file_list_sample %in% metadata$ALT_NAME)]

# Read all the processed assemblytics files
assemblytics_ALL = NULL
for (i in file_list){
  
  print(i)
  # Read individual processed assemblytics file
  assemblytics = read.table(paste0(file_path,i), stringsAsFactors = F, header = T)
  
  # Combine
  assemblytics_ALL = rbind.data.frame(assemblytics_ALL, assemblytics, stringsAsFactors = F)
  
}

#sample_list = unique(assemblytics_ALL$sample)

#assemblytics_ALL_pseudohap_combined = combinePseudohaps(sample_list, assemblytics_ALL)

# Turn off scientific notation
options(scipen=999)
write.table(assemblytics_ALL, paste0(dir, "/discovery/assemblytics_combined_results.txt"), col.names = T, row.names = F, quote = F, sep = '\t')



